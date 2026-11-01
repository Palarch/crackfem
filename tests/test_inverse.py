import numpy as np
import pytest

from crackfem.inverse import estimate_eta, eta_decompression, is_monotone_nonincreasing


def w_pred(eta):
    return 1e-3 * max(0.0, 0.6 - eta)


@pytest.mark.parametrize("eta_true", [0.15, 0.4, 0.55])
def test_round_trip(eta_true):
    est = estimate_eta(w_pred, w_pred(eta_true))
    assert est.status == "ok"
    assert est.eta == pytest.approx(eta_true, abs=1e-6)


def test_no_crack_returns_a_bound_not_a_value():
    est = estimate_eta(w_pred, 0.0)
    assert est.status == "no_crack_bound"
    assert est.eta is None
    assert est.eta_lower_bound == pytest.approx(0.6, abs=1e-5)


def test_crack_wider_than_the_model_can_explain():
    est = estimate_eta(w_pred, 1.0e-3)
    assert est.status == "exceeds_model"
    assert est.eta is None


def test_intact_girder_that_already_cracks_is_reported_as_inconsistent():
    with pytest.raises(ValueError):
        estimate_eta(lambda e: 1e-3 * (1.2 - e), 0.1e-3)


def test_decompression_eta_and_monotonicity():
    assert eta_decompression(w_pred) == pytest.approx(0.6, abs=1e-5)
    assert is_monotone_nonincreasing(w_pred, np.linspace(0.02, 1.0, 50))
    assert not is_monotone_nonincreasing(lambda e: e, np.linspace(0.02, 1.0, 50))


def counted(w):
    calls = [0]

    def g(eta):
        calls[0] += 1
        return w(eta)

    return g, calls


def w_smooth(eta):                 # nonlinear, no plateau kink in its derivative at the root
    return 1e-3 * (0.6 - eta) ** 2 if eta < 0.6 else 0.0


def w_jump(eta):                   # width jumps from 0.0668 mm to zero at decompression
    return 0.0668e-3 + 1e-3 * (0.6 - eta) if eta < 0.6 else 0.0


def test_non_finite_measured_width_is_refused():
    for bad in (float("nan"), float("inf"), -float("inf")):
        with pytest.raises(ValueError):
            estimate_eta(w_pred, bad)


@pytest.mark.parametrize("w_meas", [-1e-4, -1e-12, 1e-10, 1e-9])
def test_widths_at_or_below_a_nanometre_are_no_crack(w_meas):
    est = estimate_eta(w_pred, w_meas)
    assert est.status == "no_crack_bound"
    assert est.eta is None
    assert est.eta_lower_bound == pytest.approx(0.6, abs=1e-5)


def test_just_above_a_nanometre_is_still_a_crack():
    est = estimate_eta(w_pred, 2e-9)
    assert est.status == "ok"
    assert est.eta == pytest.approx(0.6 - 2e-6, abs=1e-7)


def test_decompression_returns_the_lower_end_when_the_lowest_eta_is_already_uncracked():
    w_low = lambda e: 1e-3 * max(0.0, 0.01 - e)       # no crack at or above eta = 0.01
    assert eta_decompression(w_low, eta_lo=0.02) == 0.02


def test_decompression_refuses_an_intact_girder_that_already_cracks():
    with pytest.raises(ValueError):
        eta_decompression(lambda e: 1e-3 * (1.2 - e))


def test_decompression_returns_the_cracked_side_of_the_edge_within_xtol():
    edge = eta_decompression(w_pred, xtol=1e-7)
    assert w_pred(edge) > 1e-9
    assert 0.6 - 1e-6 - 1e-7 <= edge <= 0.6 - 1e-6 + 1e-12      # root of w = w_tiny is 0.6 - 1e-6


def test_decompression_bisection_cost():
    g, calls = counted(w_pred)
    eta_decompression(g)
    assert calls[0] <= 30


def test_decompression_rejects_a_nan_forward_model():
    with pytest.raises(ValueError):
        eta_decompression(lambda e: float("nan"))


def test_xtol_reaches_the_decompression_search_of_the_bound():
    fine = estimate_eta(w_pred, 0.0, xtol=1e-7)
    coarse = estimate_eta(w_pred, 0.0, xtol=1e-3)
    assert coarse.evaluations < fine.evaluations
    assert coarse.eta_lower_bound == pytest.approx(0.6, abs=1e-3)


@pytest.mark.parametrize("eta_true", [0.1, 0.3, 0.5])
def test_smooth_nonlinear_round_trip_to_xtol(eta_true):
    est = estimate_eta(w_smooth, w_smooth(eta_true), xtol=1e-7)
    assert est.status == "ok"
    assert est.eta == pytest.approx(eta_true, abs=1e-7)


def test_a_jump_in_the_forward_map_gives_a_bound_not_a_point_estimate():
    est = estimate_eta(w_jump, 0.02e-3)
    assert est.status == "no_crack_bound"
    assert est.eta is None
    assert est.eta_lower_bound == pytest.approx(0.6, abs=1e-5)


def test_widths_the_jump_map_can_produce_still_round_trip():
    est = estimate_eta(w_jump, w_jump(0.3))
    assert est.status == "ok"
    assert est.eta == pytest.approx(0.3, abs=1e-6)


def test_evaluations_are_reported_for_an_ok_case():
    g, calls = counted(w_pred)
    est = estimate_eta(g, w_pred(0.4))
    assert est.status == "ok"
    assert est.evaluations == calls[0]
    assert est.evaluations == est.iterations + 3      # f(hi), f(lo), the iterates, the residual check


def test_evaluations_are_reported_for_the_bound_case():
    g, calls = counted(w_pred)
    est = estimate_eta(g, 0.0)
    assert est.evaluations == calls[0] <= 30


def test_monotone_check_needs_a_strictly_ascending_grid():
    for bad in ([0.5, 0.4], [0.1, 0.1, 0.2], [0.1, float("nan")]):
        with pytest.raises(ValueError):
            is_monotone_nonincreasing(w_pred, bad)


def test_monotone_check_tolerance_scales_with_the_width():
    grid = [0.2, 0.4, 0.6, 0.8]
    noisy = dict(zip(grid, [4e-4, 2e-4, 2e-4 + 1e-14, 0.0]))       # a 1e-14 round-off rise
    broken = dict(zip(grid, [4e-4, 2e-4, 2e-4 + 1e-6, 0.0]))       # a real rise, 1e-6 m
    assert is_monotone_nonincreasing(lambda e: noisy[e], grid)     # below 1e-9 * max|w| = 4e-13
    assert not is_monotone_nonincreasing(lambda e: noisy[e], grid, tol=0.0)
    assert not is_monotone_nonincreasing(lambda e: broken[e], grid)
    assert is_monotone_nonincreasing(lambda e: 0.0, grid)          # all zero: absolute 1e-15 floor
