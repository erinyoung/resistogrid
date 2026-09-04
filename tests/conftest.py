from pathlib import Path
import pytest


@pytest.fixture
def test_data_dir() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture
def amrfinder_plus_file(test_data_dir: Path) -> Path:
    return test_data_dir / "amrfinder_plus.tsv"


@pytest.fixture
def amrfinder_legacy_file(test_data_dir: Path) -> Path:
    return test_data_dir / "amrfinder_legacy.tsv"