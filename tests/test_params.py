import re
from pathlib import Path

from crackfem.params import SOURCES, Inputs, field_names

SOURCES_MD = Path(__file__).resolve().parents[1] / "docs" / "sources.md"


def _flatten(value):
    if isinstance(value, tuple):
        for item in value:
            yield from _flatten(item)
    else:
        yield value


def test_every_input_has_a_source_id():
    assert field_names() == set(SOURCES)


def test_every_source_id_is_a_row_of_sources_md():
    text = SOURCES_MD.read_text(encoding="utf-8")
    ids = re.findall(r"^\|\s*(S\d\d)\s*\|", text, flags=re.MULTILINE)
    assert len(ids) == len(set(ids))
    assert set(SOURCES.values()) <= set(ids)


def test_all_numeric_inputs_are_positive():
    inp = Inputs()
    for name in field_names():
        assert all(v > 0 for v in _flatten(getattr(inp, name))), name
