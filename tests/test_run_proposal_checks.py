"""run_proposal_checks: the hand calculation and the width-response rows (Expected Outcome 2)."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from crackfem.girder1d import Girder1D

_spec = importlib.util.spec_from_file_location(
    "run_proposal_checks", Path(__file__).resolve().parents[1] / "scripts" / "run_proposal_checks.py")
rpc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rpc)


def test_midspan_ratio_is_the_smaller_of_the_two_elements_at_midspan():
    assert rpc.midspan_ratio([1.0, 0.8, 0.3, 0.4, 0.9, 1.0]) == 0.3


def test_hand_calculation_matches_the_code_and_the_closed_form():
    h = rpc.hand_calc()
    Le, EI = h["Le_m"], h["EI_kNm2"]
    k = np.array(h["k_e_kN_m"])
    assert k[0][0] == pytest.approx(12 * EI / Le**3, rel=1e-12)
    assert k[1][1] == pytest.approx(4 * EI / Le, rel=1e-12)
    assert np.allclose(k, k.T)
    q, L = h["q_kN_m"], h["span_m"]
    assert h["f_e_kN"][0] == pytest.approx(-q * Le / 2, rel=1e-12)          # q acts downward
    assert sum(h["reactions_kN"]) == pytest.approx(q * L, rel=1e-12)
    assert h["v_mid_fe_mm"] == pytest.approx(h["v_mid_exact_mm"], rel=1e-9)
    assert h["node2_block"][0][1] == pytest.approx(0.0, abs=1e-9 * EI)       # the 6EI/L^2 terms cancel


def test_stage1_rows_start_intact_and_soften_as_the_crack_widens():
    g = Girder1D(n_el=10)
    rows = rpc.stage1_rows(g, [1.0, 0.4, 0.2])
    assert rows[0]["w_mm"] == 0.0
    assert rows[0]["alpha_equivalent"] == pytest.approx(1.0, rel=1e-12)
    assert rows[0]["EI_ratio_midspan"] == pytest.approx(1.0, rel=1e-9)
    assert rows[0]["w_mm"] < rows[1]["w_mm"] < rows[2]["w_mm"]
    assert rows[0]["deflection_truck_mm"] < rows[1]["deflection_truck_mm"] < rows[2]["deflection_truck_mm"]
    assert rows[2]["EI_ratio_midspan"] < rows[1]["EI_ratio_midspan"] < 1.0


def test_stage1_band_orders_each_response_from_low_to_high():
    g = Girder1D(n_el=10)
    case = {"name": "c", "w_meas_mm": 0.16, "u95_mm": 0.04, "eta_mc": {"p05": 0.34, "p50": 0.36, "p95": 0.38}}
    b = rpc.stage1_case_band(g, case, g.indicator(1.0)["deflection_truck_mm"])
    assert (b["w"], b["u"]) == (0.16, 0.04)
    for key in ("deflection_truck_mm", "alpha_equivalent", "EI_ratio_midspan"):
        assert b[key]["lo"] <= b[key]["mid"] <= b[key]["hi"], key
    assert b["deflection_truck_mm"]["mid"] == pytest.approx(g.indicator(0.36)["deflection_truck_mm"], rel=1e-12)


def test_end_points_solve_again_at_w_minus_u_and_w_plus_u():
    g = Girder1D(n_el=10)
    case = {"name": "c", "gsd_mm": 0.02, "w_meas_mm": 0.16, "u95_mm": 0.04, "eta_point": 0.36,
            "eta_mc": {"p05": 0.34, "p95": 0.38}}
    (r,) = rpc.end_point_rows(g, [case])
    assert r["status"] == ["ok", "ok"]
    assert r["eta_at_w_plus_u"] < r["eta_at_w_minus_u"]                    # a wider crack means less prestress
    assert g.w_pred(r["eta_at_w_plus_u"]) == pytest.approx(0.20e-3, rel=1e-3)
    assert g.w_pred(r["eta_at_w_minus_u"]) == pytest.approx(0.12e-3, rel=1e-3)
    assert (r["mc_p05"], r["mc_p95"]) == (0.34, 0.38)
