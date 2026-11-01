from dataclasses import replace

import numpy as np
import pytest

from crackfem.girder1d import Girder1D
from crackfem.inverse import estimate_eta, eta_decompression, is_monotone_nonincreasing
from crackfem.params import Inputs
from crackfem.solve import solve


@pytest.fixture(scope="module")
def g():
    return Girder1D()


def test_total_dead_moment_is_stage_a_plus_stage_b(g):
    L = g.inp.span
    assert g.M_dead_mid == pytest.approx((g.q_a + g.q_b) * L**2 / 8, rel=1e-9)


def test_intact_girder_shows_no_crack_at_rest(g):
    assert g.w_pred(1.0) == 0.0


def test_predicted_width_does_not_increase_with_eta(g):
    assert is_monotone_nonincreasing(g.w_pred, np.linspace(0.02, 1.0, 40))


def test_round_trip_through_the_real_section(g):
    eta_true = round(eta_decompression(g.w_pred) - 0.15, 3)
    est = estimate_eta(g.w_pred, g.w_pred(eta_true))
    assert est.status == "ok"
    assert est.eta == pytest.approx(eta_true, abs=1e-5)


def test_intact_service_response_is_the_linear_uncracked_solution(g):
    ind = g.indicator(1.0)
    lin = solve(g.model_service, g.EI_I)
    assert ind["deflection_service_mm"] == pytest.approx(-lin.D[2 * g.mid_node] * 1e3, rel=1e-9)
    assert all(r == pytest.approx(1.0) for r in ind["EI_ratio"])
    assert ind["visible_crack_zone_m"] == 0.0


def test_damaged_girder_is_softer_but_still_one_pass(g):
    eta = round(eta_decompression(g.w_pred) - 0.15, 3)
    ind, ref = g.indicator(eta), g.indicator(1.0)
    assert ind["deflection_service_mm"] > ref["deflection_service_mm"]
    assert any(ind["precracked"])
    assert min(ind["EI_ratio"]) < 1.0
    assert ind["iterations_service"] == 1
    assert 0.0 < ind["visible_crack_zone_m"] < g.inp.span


def test_camber_is_upward_and_stage_a_deflection_downward(g):
    ind = g.indicator(1.0)
    assert ind["camber_mm"] > 0.0
    assert ind["deflection_stage_a_mm"] > 0.0


@pytest.mark.parametrize("eta", [0.02, 0.06])
def test_very_low_eta_runs_and_returns_finite_deflections(g, eta):
    # at very low eta the total-strain zero of a staged section can lie outside the girder
    ind = g.indicator(eta)
    for k in ("deflection_dead_mm", "deflection_service_mm", "deflection_truck_mm"):
        assert np.isfinite(ind[k])


def test_precracked_zone_uses_the_moving_truck_envelope(g):
    assert np.all(g.M_service_env_el >= g.M_service_el - 1e-9 * g.M_service)
    assert np.any(g.M_service_env_el > 1.001 * g.M_service_el)     # off midspan the truck moves closer
    assert not any(g.precracked_zone(1.0))


def test_branson_equals_ec2_on_the_intact_girder(g):
    # at eta = 1 no element passes Mcr, so both methods keep the uncracked stiffness
    b, e = g.indicator(1.0, method="branson"), g.indicator(1.0)
    for k in ("deflection_dead_mm", "deflection_service_mm", "deflection_truck_mm"):
        assert b[k] == pytest.approx(e[k], rel=1e-12)


@pytest.mark.parametrize("eta", [0.407, 1.0])
def test_kappa_curves_start_at_the_end_of_stage_a(g, eta):
    M_A = float(g.M_A_nd[g.mid_node])
    M_grid = np.linspace(M_A, 1.3 * g.M_service, 25)
    kc = g.kappa_curves(eta, M_grid)
    mdec = g.mcr(eta, g.sec_nd[g.mid_node], 0.0)
    cracked = np.array(kc["cracked"])
    assert np.all(np.isnan(cracked[M_grid < mdec]))
    assert np.all(np.isfinite(cracked[M_grid >= mdec]))
    for k in ("uncracked", "ec2_sound", "ec2_precracked"):     # each curve measured from its own kappa(M_A)
        assert kc[k][0] == 0.0


def test_as_built_strand_count_reproduces_the_designed_girder(g):
    a = Girder1D(n_strands=g.design.n_strands)
    assert a.design == g.design
    for eta in (0.2, 0.4, 0.6):
        assert a.w_pred(eta) == g.w_pred(eta)
    assert a.indicator(0.4) == g.indicator(0.4)


def test_as_built_strand_count_is_kept_when_an_input_changes(g):
    # an existing girder does not get new strands when an input is perturbed
    inp = replace(Inputs(), girder_spacing=2.6)
    redesigned = Girder1D(inp)
    assert redesigned.design.n_strands != g.design.n_strands      # the design search would change it
    a = Girder1D(inp, n_strands=g.design.n_strands)
    assert a.design.n_strands == g.design.n_strands
    assert a.design.A_p == g.design.n_strands * inp.strand_area
    assert a.design.M_service == a.M_service                     # soffit stress recorded at this girder's M_service
    assert a.design.soffit_stress_service > 0.0                  # fewer strands than this load needs: tension
