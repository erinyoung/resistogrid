from resistogrid.io import expand_file_inputs


def test_expand_file_inputs(tmp_path):
    sub_dir = tmp_path / "batch"
    sub_dir.mkdir()

    file1 = tmp_path / "sample1.tsv"
    file2 = sub_dir / "sample2.txt"
    file_ignored = sub_dir / "sample3.pdf"

    file1.touch()
    file2.touch()
    file_ignored.touch()

    expanded = expand_file_inputs([tmp_path])
    assert file1 in expanded
    assert file2 in expanded
    assert file_ignored not in expanded