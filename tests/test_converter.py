import pytest
import pandas as pd
from resistogrid.converter import read_and_normalize, process_amrfinder_files


@pytest.mark.parametrize("fixture_name", ["amrfinder_plus_file", "amrfinder_legacy_file"])
def test_read_and_normalize_schemas(fixture_name, request):
    file_path = request.getfixturevalue(fixture_name)
    df = read_and_normalize(file_path)

    assert not df.empty
    # Verify canonical normalized column names exist regardless of input schema
    assert "sample" in df.columns
    assert "symbol" in df.columns
    assert "coverage" in df.columns
    assert "identity" in df.columns


def test_process_combined_schemas_matrix(amrfinder_plus_file, amrfinder_legacy_file, tmp_path):
    output_matrix = tmp_path / "combined_matrix.tsv"

    process_amrfinder_files(
        files=[amrfinder_plus_file, amrfinder_legacy_file],
        output=output_matrix,
        matrix_value="binary",
        group_by="hierarchical",
    )

    assert output_matrix.is_file()

    # Load produced matrix and check structure
    result_df = pd.read_csv(output_matrix, sep="\t", index_col="Sample")
    assert not result_df.empty
    assert len(result_df) >= 2  # Samples from both schema files present