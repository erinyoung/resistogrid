import json
import re
import sys
import urllib.request
from importlib.resources import files
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

NCBI_HIERARCHY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt"
)
FETCH_COLUMNS = ["node_id", "parent_node_id", "symbol", "db_version"]
CSV_COLUMNS = ["node_id", "parent_node_id", "symbol"]


def _find_family_symbol(node_id: str, node_dict: Dict[str, Dict[str, str]]) -> str:
    """Traverse the hierarchy tree upwards to locate the parent family symbol."""
    curr_id = node_id
    visited = set()

    while curr_id in node_dict and curr_id not in visited:
        visited.add(curr_id)
        row = node_dict[curr_id]

        # Priority 1: Family node with explicit _fam suffix
        if curr_id.endswith("_fam"):
            return curr_id[:-4]

        # Priority 2: Intermediate parent marked as display_parent == '1'
        if row.get("display_parent") == "1" and curr_id != node_id:
            sym = row.get("symbol")
            return sym if sym else curr_id

        parent_id = row.get("parent_node_id")
        if not parent_id or parent_id in ["ALL", "AMR", "STRESS", "VIRULENCE", "ROOT"]:
            sym = row.get("symbol")
            return sym if sym else curr_id

        curr_id = parent_id

    return node_id[:-4] if node_id.endswith("_fam") else node_id


def _build_mapping_from_df(df: pd.DataFrame) -> Dict[str, str]:
    """Parse node hierarchy dataframe into a gene -> family mapping dictionary."""
    records = df.fillna("").to_dict(orient="records")
    node_dict = {
        row["node_id"].strip(): row
        for row in records
        if row.get("node_id")
    }

    mapping: Dict[str, str] = {}
    for row in records:
        node_id = row.get("node_id", "").strip()
        symbol = row.get("symbol", "").strip()

        if not node_id:
            continue

        if not symbol:
            symbol = node_id

        family = _find_family_symbol(node_id, node_dict)

        # Don't store a gene -> itself mapping
        if symbol and family and symbol != family:
            mapping[symbol] = family

    return mapping


def update_gene_key(target_dir: Optional[Path] = None) -> Path:
    if target_dir is None:
        target_dir = Path(__file__).resolve().parent / "data"

    target_dir.mkdir(parents=True, exist_ok=True)
    csv_file = target_dir / "ncbi_hierarchy.csv"
    version_file = target_dir / "ncbi_version.txt"

    print("Downloading latest NCBI ReferenceGeneHierarchy.txt...", file=sys.stderr)

    try:
        req = urllib.request.Request(
            NCBI_HIERARCHY_URL,
            headers={"User-Agent": "resistogrid-maintainer"},
        )

        with urllib.request.urlopen(req) as response:
            df = pd.read_csv(
                response,
                sep="\t",
                usecols=FETCH_COLUMNS,
                dtype=str,
                low_memory=False,
            )

        # Separate db_version metadata
        if "db_version" in df.columns:
            versions = df["db_version"].fillna("").str.strip()
            versions = versions[versions != ""]
            if not versions.empty:
                version_file.write_text(f"{versions.iloc[0]}\n", encoding="utf-8")

        # Save lean 3-column CSV
        df = df[CSV_COLUMNS].fillna("")
        df.to_csv(csv_file, index=False)

        size_kb = csv_file.stat().st_size / 1024
        print(
            f"Successfully saved lean hierarchy ({len(df):,} rows, {size_kb:.1f} KB) to:\n"
            f"  {csv_file}",
            file=sys.stderr,
        )

        return csv_file

    except Exception as e:
        print(
            f"Error: Failed to download or process NCBI gene hierarchy: {e}",
            file=sys.stderr,
        )
        raise


def load_gene_key(key_file: Optional[Path] = None) -> Dict[str, str]:
    """
    Load the gene hierarchy CSV and construct the gene -> family mapping.
    """
    path_to_load = None

    if key_file is not None:
        if key_file.is_file():
            path_to_load = key_file
        else:
            print(
                f"Warning: Specified key file '{key_file}' not found. "
                "Falling back to default package data.",
                file=sys.stderr,
            )

    if path_to_load is None:
        try:
            pkg_data = files("resistogrid.data").joinpath("ncbi_hierarchy.csv")
            if pkg_data.is_file():
                path_to_load = pkg_data
        except Exception:
            pass

    if path_to_load is None:
        local_default = (
            Path(__file__).resolve().parent / "data" / "ncbi_hierarchy.csv"
        )
        if local_default.is_file():
            path_to_load = local_default

    if path_to_load is not None:
        try:
            # Backward compatibility for legacy JSON formats
            if str(path_to_load).endswith(".json"):
                with open(path_to_load, "r", encoding="utf-8") as f:
                    data = json.load(f)
                return data.get("gene_to_family", data)

            df = pd.read_csv(path_to_load, dtype=str)
            return _build_mapping_from_df(df)
        except Exception as e:
            print(
                f"Warning: Could not load gene family key from '{path_to_load}': {e}",
                file=sys.stderr,
            )

    print(
        "Warning: Could not load gene family key. "
        "Using Sequence Name / regex fallback logic.",
        file=sys.stderr,
    )

    return {}


def resolve_family(
    gene_symbol: str,
    seq_name: str = "",
    key_map: Optional[Dict[str, str]] = None,
) -> str:
    """
    Return mapped family symbol for a gene/allele.
    1. Checks the key_map dictionary (built from NCBI ReferenceGeneHierarchy).
    2. Parses "X family" from Sequence Name / Element Name field.
    3. Falls back to regex pattern stripping (e.g., blaCTX-M-27 -> blaCTX-M).
    """
    if key_map is None:
        key_map = {}

    gene_symbol = str(gene_symbol).strip() if gene_symbol else ""
    seq_name = str(seq_name).strip() if seq_name and str(seq_name) != "nan" else ""

    if gene_symbol in key_map:
        return key_map[gene_symbol]

    if seq_name and (
        gene_symbol.lower().startswith("bla")
        or "beta-lactamase" in seq_name.lower()
    ):
        m = re.search(
            r"([A-Za-z0-9_\-\(\)\'\"]+)\s+family", seq_name, re.IGNORECASE
        )
        if m:
            fam_name = m.group(1).strip()
            if gene_symbol.lower().startswith("bla"):
                if not fam_name.lower().startswith("bla"):
                    fam_name = "bla" + fam_name
                elif fam_name.lower().startswith("bla") and fam_name[:3] != "bla":
                    fam_name = "bla" + fam_name[3:]
            return fam_name

    heuristic = re.sub(r"[-_]\d+.*$", "", gene_symbol)
    return heuristic if heuristic else gene_symbol