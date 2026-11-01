"""Stage-2 inputs: values and source rows (S22-S27)."""
import re
from pathlib import Path

import pytest

from crackfem.params import FT, SOURCES, Inputs

SOURCES_MD = Path(__file__).resolve().parents[1] / "docs" / "sources.md"
STAGE2 = {"overhang_width_ft": (3.0, "S25"), "barrier_width_ft": (1.0, "S25"),
          "wheel_spacing_ft": (6.0, "S26"), "truck_clearance_ft": (10.0, "S26"),
          "wheel_curb_clearance_ft": (2.0, "S26"), "lane_width_ft": (12.0, "S26"),
          "poisson_concrete": (0.2, "S24"), "diaphragm_depth": (1.10, "S27"),
          "diaphragm_thickness": (0.25, "S27")}


@pytest.mark.parametrize("name", sorted(STAGE2))
def test_stage2_input_value_and_source(name):
    value, sid = STAGE2[name]
    assert getattr(Inputs(), name) == value
    assert SOURCES[name] == sid


def test_new_source_rows_exist_and_s22_is_filled():
    text = SOURCES_MD.read_text(encoding="utf-8")
    ids = re.findall(r"^\|\s*(S\d\d)\s*\|", text, flags=re.MULTILINE)
    assert {"S22", "S23", "S24", "S25", "S26", "S27"} <= set(ids)
    s22 = next(line for line in text.splitlines() if line.startswith("| S22 |"))
    assert "Example 5.5" in s22 and "local PDF" not in s22


def test_truck_fits_its_clearance_width():
    inp = Inputs()
    assert 2 * inp.wheel_curb_clearance_ft + inp.wheel_spacing_ft == pytest.approx(inp.truck_clearance_ft)
    assert inp.overhang_width_ft * FT > inp.barrier_width_ft * FT
