import sys
from pathlib import Path
from typing import List
import pandas as pd


def expand_file_inputs(paths: List[Path]) -> List[Path]:
    """
    Expand input paths to collect all valid files.
    If a path is a directory, recursively collects all .tsv and .txt files within it.
    """
    resolved_files: List[Path] = []
    for path in paths:
        if path.is_dir():
            dir_files = [
                p
                for p in path.rglob("*")
                if p.is_file() and p.suffix.lower() in [".tsv", ".txt", ".csv"]
            ]
            resolved_files.extend(dir_files)
        elif path.is_file():
            resolved_files.append(path)
        else:
            print(f"Warning: Input path '{path}' does not exist.", file=sys.stderr)

    # Return unique, sorted path list
    return sorted(list(set(resolved_files)))


def read_tsv_raw(file_path: Path) -> pd.DataFrame:
    """Safely load a raw tab-delimited AMRFinder file into a DataFrame as strings."""
    try:
        return pd.read_csv(file_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"Warning: Failed to read file '{file_path}': {e}", file=sys.stderr)
        return pd.DataFrame()


def export_matrix(matrix: pd.DataFrame, output_path: Path) -> None:
    """Write the processed matrix to disk. Uses comma separation for .csv and tab for all others."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sep = "," if output_path.suffix.lower() == ".csv" else "\t"
    matrix.to_csv(output_path, sep=sep)