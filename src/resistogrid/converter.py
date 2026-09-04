import sys
from pathlib import Path
from typing import Dict, List, Optional, Set
import pandas as pd

from resistogrid.key import load_gene_key, resolve_family

SCHEMA_MAPPINGS: Dict[str, Dict[str, str]] = {
    "plus": {
        "sample": "Name",
        "symbol": "Gene symbol",
        "name": "Sequence name",
        "type": "Element type",
        "method": "Method",
        "target_len": "Target length",
        "ref_len": "Reference sequence length",
        "coverage": "% Coverage of reference sequence",
        "identity": "% Identity to reference sequence",
        "aln_len": "Alignment length",
    },
    "legacy": {
        "sample": "Name",
        "symbol": "Element symbol",
        "name": "Element name",
        "type": "Type",
        "method": "Method",
        "target_len": "Target length",
        "ref_len": "Reference sequence length",
        "coverage": "% Coverage of reference",
        "identity": "% Identity to reference",
        "aln_len": "Alignment length",
    },
}


def detect_schema(columns: List[str]) -> Dict[str, str]:
    cols_set = set(columns)
    if "Gene symbol" in cols_set:
        return SCHEMA_MAPPINGS["plus"]
    elif "Element symbol" in cols_set:
        return SCHEMA_MAPPINGS["legacy"]
    else:
        raise ValueError("Unrecognized file header schema.")


def read_and_normalize(file_path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"Warning: Could not read '{file_path}': {e}", file=sys.stderr)
        return pd.DataFrame()

    if df.empty:
        return pd.DataFrame()

    schema = detect_schema(df.columns.tolist())
    rename_dict = {v: k for k, v in schema.items() if v in df.columns}
    df = df.rename(columns=rename_dict)

    if "sample" not in df.columns or df["sample"].dropna().empty:
        df["sample"] = file_path.stem

    numeric_cols = ["target_len", "ref_len", "coverage", "identity", "aln_len"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def apply_filters(
    df: pd.DataFrame,
    gene_filter: Optional[str] = None,
    element_types: Optional[List[str]] = None,
    methods: Optional[List[str]] = None,
    min_target_len: Optional[int] = None,
    min_ref_len: Optional[int] = None,
    min_cov: Optional[float] = None,
    min_ident: Optional[float] = None,
    min_aln_len: Optional[int] = None,
) -> pd.DataFrame:
    if df.empty:
        return df

    if gene_filter:
        pattern = gene_filter.lower()
        symbol_match = df["symbol"].fillna("").str.lower().str.contains(pattern)
        name_match = (
            df["name"].fillna("").str.lower().str.contains(pattern)
            if "name" in df.columns
            else False
        )
        df = df[symbol_match | name_match]

    if element_types and "type" in df.columns:
        types_upper = {t.upper() for t in element_types}
        df = df[df["type"].fillna("").str.upper().isin(types_upper)]

    if methods and "method" in df.columns:
        methods_upper = {m.upper() for m in methods}
        df = df[df["method"].fillna("").str.upper().isin(methods_upper)]

    if min_target_len is not None and "target_len" in df.columns:
        df = df[df["target_len"] >= min_target_len]
    if min_ref_len is not None and "ref_len" in df.columns:
        df = df[df["ref_len"] >= min_ref_len]
    if min_cov is not None and "coverage" in df.columns:
        df = df[df["coverage"] >= min_cov]
    if min_ident is not None and "identity" in df.columns:
        df = df[df["identity"] >= min_ident]
    if min_aln_len is not None and "aln_len" in df.columns:
        df = df[df["aln_len"] >= min_aln_len]

    return df


def pivot_matrix(
    df: pd.DataFrame,
    all_samples: Set[str],
    matrix_value: str,
    group_by: str,
    key_map: Dict[str, str],
) -> pd.DataFrame:
    if df.empty:
        res = pd.DataFrame(index=sorted(all_samples))
        res.index.name = "Sample"
        return res

    df = df.copy()

    # Pass symbol and optional sequence name to family resolver
    def _get_family(row: pd.Series) -> str:
        symbol = str(row.get("symbol", ""))
        name = str(row.get("name", "")) if "name" in row else ""
        return resolve_family(symbol, seq_name=name, key_map=key_map)

    df["family"] = df.apply(_get_family, axis=1)

    if group_by == "hierarchical":
        all_fams = sorted(df["family"].unique())
        fam_to_genes = {
            fam: sorted(df[df["family"] == fam]["symbol"].unique())
            for fam in all_fams
        }

        rows = []
        for sample in sorted(all_samples):
            sample_df = df[df["sample"] == sample]
            row = {"Sample": sample}

            for fam in all_fams:
                fam_df = sample_df[sample_df["family"] == fam]
                matched_genes = sorted(fam_df["symbol"].unique())
                has_fam = len(matched_genes) > 0

                # Family aggregated presence/metric column
                if matrix_value == "detected":
                    row[f"{fam}_fam"] = "Detected" if has_fam else "Not Detected"
                elif matrix_value == "binary":
                    row[f"{fam}_fam"] = 1 if has_fam else 0
                elif matrix_value == "count":
                    row[f"{fam}_fam"] = len(fam_df)
                elif matrix_value in ["identity", "coverage"]:
                    row[f"{fam}_fam"] = fam_df[matrix_value].max() if has_fam else 0.0

                # Comma-separated list of detected allele symbols in this family
                row[f"{fam}_fam_elements"] = ", ".join(matched_genes)

                # Individual gene columns
                for gene in fam_to_genes[fam]:
                    gene_df = fam_df[fam_df["symbol"] == gene]
                    has_gene = len(gene_df) > 0

                    if matrix_value == "detected":
                        row[gene] = "Detected" if has_gene else "Not Detected"
                    elif matrix_value == "binary":
                        row[gene] = 1 if has_gene else 0
                    elif matrix_value == "count":
                        row[gene] = len(gene_df)
                    elif matrix_value in ["identity", "coverage"]:
                        row[gene] = gene_df[matrix_value].max() if has_gene else 0.0

            rows.append(row)

        res_df = pd.DataFrame(rows).set_index("Sample")
        return res_df

    df["family_symbol"] = df["family"] + "___" + df["symbol"]
    target_col = {
        "symbol": "symbol",
        "family": "family",
        "family_symbol": "family_symbol",
    }.get(group_by, "symbol")

    if matrix_value == "detected":
        matrix = df.pivot_table(
            index="sample",
            columns=target_col,
            values="identity",
            aggfunc=lambda x: "Detected" if len(x) > 0 else "Not Detected",
            fill_value="Not Detected",
        )
    elif matrix_value == "binary":
        matrix = df.pivot_table(
            index="sample",
            columns=target_col,
            values="identity",
            aggfunc=lambda x: 1 if len(x) > 0 else 0,
            fill_value=0,
        )
    elif matrix_value == "count":
        matrix = df.pivot_table(
            index="sample",
            columns=target_col,
            values="identity",
            aggfunc="count",
            fill_value=0,
        )
    elif matrix_value in ["identity", "coverage"]:
        val_col = "identity" if matrix_value == "identity" else "coverage"
        matrix = df.pivot_table(
            index="sample",
            columns=target_col,
            values=val_col,
            aggfunc="max",
            fill_value=0,
        )
    else:
        raise ValueError(f"Unsupported matrix value mode: {matrix_value}")

    fill_val = "Not Detected" if matrix_value == "detected" else 0
    matrix = matrix.reindex(sorted(all_samples), fill_value=fill_val)
    matrix.index.name = "Sample"
    return matrix


def process_amrfinder_files(
    files: List[Path],
    output: Path,
    matrix_value: str = "binary",
    group_by: str = "hierarchical",
    sankey_out: Optional[Path] = None,
    key_file: Optional[Path] = None,
    gene_filter: Optional[str] = None,
    element_types: Optional[List[str]] = None,
    methods: Optional[List[str]] = None,
    min_target_len: Optional[int] = None,
    min_ref_len: Optional[int] = None,
    min_cov: Optional[float] = None,
    min_ident: Optional[float] = None,
    min_aln_len: Optional[int] = None,
) -> None:
    key_map = load_gene_key(key_file=key_file)
    normalized_dfs: List[pd.DataFrame] = []
    all_samples: Set[str] = set()

    for file_path in files:
        df = read_and_normalize(file_path)
        if not df.empty:
            all_samples.update(df["sample"].dropna().unique())
            filtered_df = apply_filters(
                df=df,
                gene_filter=gene_filter,
                element_types=element_types,
                methods=methods,
                min_target_len=min_target_len,
                min_ref_len=min_ref_len,
                min_cov=min_cov,
                min_ident=min_ident,
                min_aln_len=min_aln_len,
            )
            normalized_dfs.append(filtered_df)
        else:
            all_samples.add(file_path.stem)

    combined_df = (
        pd.concat(normalized_dfs, ignore_index=True)
        if normalized_dfs
        else pd.DataFrame()
    )

    if sankey_out and not combined_df.empty:
        sankey_df = combined_df[["sample", "symbol"]].drop_duplicates().copy()
        if "name" in combined_df.columns:
            sankey_df["name"] = combined_df["name"]

        sankey_df["family"] = sankey_df.apply(
            lambda r: resolve_family(
                r["symbol"], seq_name=r.get("name", ""), key_map=key_map
            ),
            axis=1,
        )
        sankey_df = sankey_df[["sample", "family", "symbol"]].rename(
            columns={"sample": "Source", "family": "Family", "symbol": "Target"}
        )
        sankey_out.parent.mkdir(parents=True, exist_ok=True)
        sankey_df.to_csv(sankey_out, sep="\t", index=False)
        print(f"Sankey long table written to: {sankey_out}")

    matrix = pivot_matrix(combined_df, all_samples, matrix_value, group_by, key_map)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output, sep="\t")