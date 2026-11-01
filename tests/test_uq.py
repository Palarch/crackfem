import inspect
import math

import numpy as np
import pytest
from scipy.stats import truncnorm

from crackfem.params import Inputs
from crackfem.uq import (DEFAULT_SEED, below_zero_fraction, monte_carlo_eta, sample_widths,
                         summarize, tornado)


def w_pred(eta):
    return 1e-3 * max(0.0, 0.6 - eta)


def test_samples_are_positive_reproducible_and_sized():
    a = sample_widths(0.2e-3, 0.1e-3, 500, seed=1)
    b = sample_widths(0.2e-3, 0.1e-3, 500, seed=1)
    assert a.size == 500 and (a > 0).all()
    np.testing.assert_array_equal(a, b)


def test_default_seed_is_shared_by_both_entry_points():
    assert inspect.signature(sample_widths).parameters["seed"].default == DEFAULT_SEED
    assert inspect.signature(monte_carlo_eta).parameters["seed"].default == DEFAULT_SEED


def test_truncation_matches_the_truncated_normal():
    w, u95, n = 0.1e-3, 0.2e-3, 20000          # sigma = 0.1e-3, so w / sigma = 1
    sigma = u95 / 2.0
    x = sample_widths(w, u95, n, seed=11)
    assert (x > 0).all()
    dist = truncnorm(-w / sigma, np.inf, loc=w, scale=sigma)
    se_mean = dist.std() / math.sqrt(n)
    assert abs(x.mean() - dist.mean()) < 3 * se_mean
    q = dist.ppf(0.05)
    se_q = math.sqrt(0.05 * 0.95 / n) / dist.pdf(q)
    assert abs(np.percentile(x, 5) - q) < 3 * se_q


def test_below_zero_fraction_is_the_normal_mass_under_zero():
    assert below_zero_fraction(0.1e-3, 0.2e-3) == pytest.approx(0.158655, abs=1e-5)   # w/sigma = 1
    assert below_zero_fraction(0.2e-3, 0.1e-3) == pytest.approx(3.167e-5, rel=1e-3)   # w/sigma = 4


@pytest.mark.parametrize("w, u", [(float("nan"), 0.1e-3), (0.2e-3, float("nan")),
                                  (float("inf"), 0.1e-3), (0.2e-3, float("inf")),
                                  (0.0, 0.1e-3), (-0.1e-3, 0.1e-3),
                                  (0.2e-3, 0.0), (0.2e-3, -0.1e-3)])
def test_sample_widths_rejects_invalid_inputs(w, u):
    with pytest.raises(ValueError):
        sample_widths(w, u, 10)


@pytest.mark.parametrize("w, u", [(float("nan"), 0.1e-3), (0.2e-3, float("nan")),
                                  (0.0, 0.1e-3), (0.2e-3, -0.1e-3)])
def test_monte_carlo_rejects_invalid_inputs(w, u):
    with pytest.raises(ValueError):
        monte_carlo_eta(w_pred, w, u, n=10)


def test_zero_uncertainty_collapses_to_the_point_estimate():
    mc = monte_carlo_eta(w_pred, 0.2e-3, 0.0, n=20)
    assert np.allclose(mc.etas, 0.4, atol=1e-6)
    assert mc.status.size == mc.widths.size == 20 and set(mc.status) == {"ok"}
    assert mc.status_counts == {"ok": 20} and mc.below_zero_fraction == 0.0


def test_zero_uncertainty_solves_once():
    calls = []

    def counted(eta):
        calls.append(eta)
        return w_pred(eta)

    monte_carlo_eta(counted, 0.2e-3, 0.0, n=500)
    one = len(calls)
    calls.clear()
    monte_carlo_eta(counted, 0.2e-3, 0.0, n=5)
    assert len(calls) == one


def test_linear_map_gives_the_analytic_band():
    # eta = 0.6 - w/1e-3, sigma_w = 0.05e-3, so p95 - p05 = 2 * 1.6449 * 0.05 = 0.16449
    mc = monte_carlo_eta(w_pred, 0.2e-3, 0.1e-3, n=20000, seed=7)
    s = mc.summary()
    assert s["p50"] == pytest.approx(0.4, abs=0.002)
    assert s["p95"] - s["p05"] == pytest.approx(0.16449, rel=0.02)
    assert mc.status_counts == {"ok": 20000}
    assert s["n_total"] == 20000 and s["n_dropped"] == 0 and s["censored"] is False


def test_dropped_statuses_are_reported_and_flag_censoring():
    n = 400
    # w_meas near the model's widest width: about half the draws exceed it -> "exceeds_model"
    mc = monte_carlo_eta(w_pred, 0.58e-3, 0.1e-3, n=n)
    assert "exceeds_model" in mc.status_counts
    assert sum(mc.status_counts.values()) == n and mc.status.size == n
    s = mc.summary()
    assert s["n_total"] == n and s["n_dropped"] > 0
    assert s["n_dropped"] == mc.status_counts["exceeds_model"]
    assert s["n"] == mc.etas.size == n - s["n_dropped"]
    assert s["censored"] is True and s["status_counts"] == mc.status_counts


def test_summary_is_none_when_no_sample_is_ok():
    mc = monte_carlo_eta(w_pred, 0.7e-3, 0.0, n=10)       # wider than the model can produce
    assert mc.etas.size == 0 and mc.status_counts == {"exceeds_model": 10}
    assert mc.summary() is None


def test_summary_percentiles_are_ordered():
    s = summarize(np.random.default_rng(3).normal(size=1000))
    assert s["p05"] < s["p50"] < s["p95"]


def test_tornado_sorts_by_swing():
    # the larger swing (slab_thickness) is alphabetically later and listed last
    metric = lambda inp: inp.girder_spacing * 1.0 + inp.slab_thickness * 100.0
    rows = tornado(metric, Inputs(), {"girder_spacing": (2.0, 2.4), "slab_thickness": (0.18, 0.22)})
    assert [r["name"] for r in rows] == ["slab_thickness", "girder_spacing"]
    assert rows[0]["swing"] == pytest.approx(4.0) and rows[1]["swing"] == pytest.approx(0.4)


@pytest.mark.parametrize("bad", ["low", "high"])
def test_tornado_marks_a_row_with_a_non_finite_metric_incomplete_and_sorts_it_last(bad):
    # NaN (e.g. a failed inverse) must not leak into max/min: the swing is NaN, the row is flagged
    # and listed after every finite row, whatever its position in the input.
    nan = float("nan")
    def metric(inp):
        if bad == "low" and inp.girder_spacing == 1.8:
            return nan
        if bad == "high" and inp.girder_spacing == 2.6:
            return nan
        return inp.girder_spacing * 1.0 + inp.slab_thickness * 100.0
    rows = tornado(metric, Inputs(), {"girder_spacing": (1.8, 2.6), "slab_thickness": (0.18, 0.22),
                                      "strand_cover": (0.035, 0.055)})
    assert [r["name"] for r in rows] == ["slab_thickness", "strand_cover", "girder_spacing"]
    assert rows[0]["swing"] == pytest.approx(4.0) and rows[1]["swing"] == pytest.approx(0.0)
    assert not rows[0].get("incomplete") and not rows[1].get("incomplete")
    assert math.isnan(rows[2]["swing"]) and rows[2]["incomplete"] is True


def test_tornado_with_a_non_finite_base_marks_every_row_incomplete():
    rows = tornado(lambda inp: float("nan"), Inputs(), {"girder_spacing": (1.8, 2.6), "slab_thickness": (0.18, 0.22)})
    assert [r["name"] for r in rows] == ["girder_spacing", "slab_thickness"]
    assert all(r["incomplete"] and math.isnan(r["swing"]) for r in rows)


def test_tornado_swing_includes_the_base_value():
    # metric peaks at the base slab thickness: low and high are equal, the swing is not zero
    metric = lambda inp: -(inp.slab_thickness - Inputs().slab_thickness) ** 2
    (row,) = tornado(metric, Inputs(), {"slab_thickness": (0.18, 0.22)})
    assert row["metric_low"] == pytest.approx(row["metric_high"])
    assert row["swing"] == pytest.approx(0.02 ** 2)
