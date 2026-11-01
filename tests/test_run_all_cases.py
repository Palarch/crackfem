"""run_all.run_case: a refusal bounds eta only when every crack is below the resolution gate."""
import importlib.util
from pathlib import Path

import pytest

from crackfem.girder1d import Girder1D

_spec = importlib.util.spec_from_file_location(
    "run_all", Path(__file__).resolve().parents[1] / "scripts" / "run_all.py")
run_all = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_all)


def _crack(i, gate):
    return {"index": i, "points": 390, "gate": gate, "max_width_px": None, "max_width_mm": None,
            "uncertainty_mm": None, "rating": None, "rating_label": None}


def _rec(gates, gsd_mm=0.02):
    return {"case": {"name": "t", "gsd_mm": gsd_mm, "eta_true": 0.36, "w_true_mm": 0.16, "width_px": 8},
            "parsed": {"gsd_mm": gsd_mm, "min_measurable_mm": 3 * gsd_mm,
                       "cracks": [_crack(i + 1, g) for i, g in enumerate(gates)]}}


@pytest.fixture(scope="module")
def g():
    return Girder1D()


def test_all_below_resolution_gives_a_refusal_with_both_bounds(g):
    out = run_all.run_case(g, _rec(["below-measurable-resolution"] * 2))
    assert out["status"] == "refused" and out["reason"] == ["below-measurable-resolution"]
    assert 0.0 < out["eta_lower_bound_95"] < out["eta_lower_bound_from_gate"] < 1.0


def test_an_unresolved_crack_gives_no_bound(g):
    out = run_all.run_case(g, _rec(["below-measurable-resolution", "no-valid-points"]))
    assert out["status"] == "unresolved"
    assert out["reason"] == ["below-measurable-resolution", "no-valid-points"]
    assert out["eta_lower_bound_from_gate"] is None and out["eta_lower_bound_95"] is None
    assert out["eta_point"] is None and out["eta_mc"] is None


def test_a_report_without_cracks_gives_no_bound(g):
    out = run_all.run_case(g, _rec([]))
    assert out["status"] == "unresolved" and out["reason"] == []
    assert out["eta_lower_bound_from_gate"] is None and out["eta_lower_bound_95"] is None


def test_stage1_inputs_leave_out_the_stage2_fields():
    """A rerun of run_all writes the same input fields as the committed stage1.json."""
    import json
    committed = json.loads((Path(__file__).resolve().parents[1] / "results" / "stage1.json")
                           .read_text(encoding="utf-8"))["inputs"]
    ins = run_all.stage1_inputs()
    assert set(ins) == set(committed)
    assert not set(run_all.STAGE2_INPUTS) & set(ins)
