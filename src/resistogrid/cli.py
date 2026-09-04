import argparse
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import List, Optional

from resistogrid.converter import process_amrfinder_files
from resistogrid.key import update_gene_key


def get_version() -> str:
    """Fetch package version from installed metadata or pyproject.toml fallback."""
    try:
        return version("resistogrid")
    except PackageNotFoundError:
        pass

    pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
    if pyproject_path.is_file():
        if sys.version_info >= (3, 11):
            import tomllib
        else:
            try:
                import tomli as tomllib
            except ImportError:
                return "0.1.0-dev"

        try:
            with open(pyproject_path, "rb") as f:
                data = tomllib.load(f)
                return data.get("project", {}).get("version", "0.1.0-dev")
        except Exception:
            pass

    return "0.1.0-dev"


__version__ = get_version()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="resistogrid",
        description="Pivot NCBI AMRFinder/AMRFinderPlus long-format outputs into sample x gene matrices.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Package Version
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program's version number and exit.",
    )

    # Database Update Command
    parser.add_argument(
        "--update-db",
        nargs="?",
        const=Path("DEFAULT"),
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "Download and optimize the latest NCBI gene family key. "
            "If DIR is omitted, saves to package data (src/resistogrid/data). "
            "If DIR is provided, saves gene_family_key.json to that directory."
        ),
    )

    # Input / Output Options
    io_group = parser.add_argument_group("Input/Output Options")
    io_group.add_argument(
        "-f",
        "--files",
        nargs="+",
        type=Path,
        required=False,
        help="One or more space-separated AMRFinder/AMRFinderPlus TSV output files.",
    )
    io_group.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("resistogrid_matrix.tsv"),
        help="Path for the output pivoted matrix TSV file.",
    )
    io_group.add_argument(
        "--sankey-out",
        type=Path,
        default=None,
        help="Optional path to output a 3-column TSV (Sample, Family, Gene) for Sankey diagrams.",
    )
    io_group.add_argument(
        "--key-file",
        type=Path,
        default=None,
        help="Optional path to a custom gene_family_key.json file instead of package defaults.",
    )

    # Matrix Value & Hierarchy Options
    matrix_group = parser.add_argument_group("Matrix Value & Hierarchy Options")
    matrix_group.add_argument(
        "--matrix-value",
        choices=["detected", "binary", "identity", "coverage", "count"],
        default="detected",
        help="Value type to populate in matrix cells: 'detected' ('Detected'/'Not Detected'), 'binary' (1/0), '%% identity', '%% coverage', or hit 'count'.",
    )
    matrix_group.add_argument(
        "--group-by",
        choices=["hierarchical", "symbol", "family", "family_symbol"],
        default="hierarchical",
        help="Matrix column grouping level: 'hierarchical' ({fam}_fam, {fam}_fam_elements, {gene}), exact 'symbol', parent 'family', or 'family_symbol' (e.g. blaOXA___blaOXA-23).",
    )

    # Filtering Options
    filter_group = parser.add_argument_group("Filtering Options")
    filter_group.add_argument(
        "-g",
        "--gene-filter",
        type=str,
        default=None,
        help="Filter results by matching substring in Gene/Element Symbol or Name (case-insensitive).",
    )
    filter_group.add_argument(
        "-t",
        "--element-type",
        nargs="+",
        type=str,
        default=None,
        help="Filter by Element Type (e.g., AMR, VIRULENCE, STRESS, POINT).",
    )
    filter_group.add_argument(
        "-m",
        "--method",
        nargs="+",
        type=str,
        default=None,
        help="Filter by identification method (e.g., BLASTX, EXACTX, PARTIALX, ALLELEX, POINTX).",
    )
    filter_group.add_argument(
        "--min-target-len",
        type=int,
        default=None,
        help="Minimum required Target length.",
    )
    filter_group.add_argument(
        "--min-ref-len",
        type=int,
        default=None,
        help="Minimum required Reference sequence length.",
    )
    filter_group.add_argument(
        "--min-cov",
        type=float,
        default=None,
        help="Minimum required %% Coverage of reference sequence (0.0 - 100.0).",
    )
    filter_group.add_argument(
        "--min-ident",
        type=float,
        default=None,
        help="Minimum required %% Identity to reference sequence (0.0 - 100.0).",
    )
    filter_group.add_argument(
        "--min-aln-len",
        type=int,
        default=None,
        help="Minimum required Alignment length.",
    )

    return parser


def main(args: Optional[List[str]] = None) -> None:
    parser = build_parser()
    parsed_args = parser.parse_args(args)

    if parsed_args.update_db is not None:
        target_dir = (
            None
            if parsed_args.update_db == Path("DEFAULT")
            else parsed_args.update_db
        )
        update_gene_key(target_dir=target_dir)
        sys.exit(0)

    if not parsed_args.files:
        parser.error(
            "the following arguments are required: -f/--files (unless running --update-db)"
        )

    missing_files = [f for f in parsed_args.files if not f.is_file()]
    if missing_files:
        print(
            f"Error: The following input file(s) do not exist:\n  "
            + "\n  ".join(str(f) for f in missing_files),
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        process_amrfinder_files(
            files=parsed_args.files,
            output=parsed_args.output,
            matrix_value=parsed_args.matrix_value,
            group_by=parsed_args.group_by,
            sankey_out=parsed_args.sankey_out,
            key_file=parsed_args.key_file,
            gene_filter=parsed_args.gene_filter,
            element_types=parsed_args.element_type,
            methods=parsed_args.method,
            min_target_len=parsed_args.min_target_len,
            min_ref_len=parsed_args.min_ref_len,
            min_cov=parsed_args.min_cov,
            min_ident=parsed_args.min_ident,
            min_aln_len=parsed_args.min_aln_len,
        )
        print(f"Success! Pivoted matrix written to: {parsed_args.output}")
    except Exception as e:
        print(f"Execution failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()