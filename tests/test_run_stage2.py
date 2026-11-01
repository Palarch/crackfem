"""run_stage2.run_case: measured, refused and unresolved cases on a tabulated forward map."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from crackfem.wtable import build_table

_spec = importlib.util.spec_from_file_location(
    "run_stage2", Path(__file__).resolve().parents[1] / "scripts" / "run_stage2.py")
run_stage2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_stage2)


def w_exact(eta):
    return 1e-3 * max(0.0, 0.6 - eta)


@pytest.fixture(scope="module")
def table():
    return build_table(w_exact, np.linspace(0.02, 1.0, 99))


def _rec(cracks, gsd_mm=0.02):
    return {"case": {"name": "t", "gsd_mm": gsd_mm, "eta_true": 0.36, "w_true_mm": 0.24, "width_px": 12},
            "parsed": {"gsd_mm": gsd_mm, "min_measurable_mm": 3 * gsd_mm, "cracks": cracks}}


def _crack(gate, w_mm=None, u_mm=None):
    return {"index": 1, "points": 390, "gate": gate, "max_width_px": None if w_mm is None else w_mm / 0.02,
            "max_width_mm": w_mm, "uncertainty_mm": u_mm, "rating": 5 if w_mm else None,
            "rating_label": "x" if w_mm else None}


def test_measured_case_band_contains_the_truth_and_the_exact_check_agrees(table):
    out = run_stage2.run_case(table, w_exact, _rec([_crack("ok", 0.24, 0.04)]))
    assert out["status"] == "measured" and out["eta_status"] == "ok"
    assert out["eta_point"] == pytest.approx(0.36, abs=1e-6)
    assert out["eta_mc"]["p05"] < 0.36 < out["eta_mc"]["p95"]
    assert out["eta_mc_exact_max_abs_diff"] < 1e-4


def test_refusal_gives_the_gate_bound_and_unresolved_gives_none(table):
    ref = run_stage2.run_case(table, w_exact, _rec([_crack("below-measurable-resolution")]))
    assert ref["status"] == "refused" and ref["eta_lower_bound_from_gate"] == pytest.approx(0.6 - 0.06, abs=1e-6)
    unr = run_stage2.run_case(table, w_exact, _rec([_crack("no-valid-points")]))
    assert unr["status"] == "unresolved" and unr["eta_lower_bound_from_gate"] is None


def test_the_generating_eta_is_labelled_as_a_1d_model_value_not_a_truth(table):
    out = run_stage2.run_case(table, w_exact, _rec([_crack("ok", 0.24, 0.04)]))
    assert out["eta_gen_1d"] == 0.36
    assert not any("true" in k for k in out)
