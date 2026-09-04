import pytest
from resistogrid.key import load_gene_key, resolve_family


def test_load_bundled_gene_key():
    key_map = load_gene_key()
    assert isinstance(key_map, dict)
    # Verify bundled database returns valid mapping
    assert len(key_map) > 0


def test_resolve_family_from_map():
    key_map = {"blaOXA-181": "blaOXA-48"}
    assert resolve_family("blaOXA-181", key_map=key_map) == "blaOXA-48"


def test_resolve_family_from_seq_name():
    # Fallback when symbol isn't in key_map but sequence name has family string
    res = resolve_family(
        gene_symbol="blaABC-1",
        seq_name="OXA-48 family class D beta-lactamase",
        key_map={},
    )
    assert res == "blaOXA-48"


def test_resolve_family_regex_heuristic():
    # Regex fallback stripping allele suffix (e.g., -27)
    assert resolve_family("blaCTX-M-27", key_map={}) == "blaCTX-M"
    assert resolve_family("tet(M)", key_map={}) == "tet(M)"