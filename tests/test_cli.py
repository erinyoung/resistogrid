import pytest
from resistogrid.cli import build_parser


def test_cli_parser_defaults():
    parser = build_parser()
    args = parser.parse_args(["-f", "input1.tsv", "input2.tsv"])

    assert args.files == [pytest.importorskip("pathlib").Path("input1.tsv"), pytest.importorskip("pathlib").Path("input2.tsv")]
    assert args.matrix_value == "detected"
    assert args.group_by == "hierarchical"


def test_cli_parser_custom_options():
    parser = build_parser()
    args = parser.parse_args(
        [
            "-f",
            "input.tsv",
            "-o",
            "out.tsv",
            "--matrix-value",
            "binary",
            "--min-cov",
            "90.0",
        ]
    )

    assert args.matrix_value == "binary"
    assert args.min_cov == 90.0