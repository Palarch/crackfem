"""Tabulated w(eta) for the Monte Carlo: exact at the knots, zero on the plateau, monotone,
and the same Monte Carlo percentiles as the exact map."""
import numpy as np
import pytest

from crackfem.uq import monte_carlo_eta
from crackfem.wtable import build_table, max_rel_error


def w_exact(eta):
    return 1e-3 * max(0.0, 0.6 - eta) ** 1.5


def test_table_is_exact_at_knots_and_zero_on_the_plateau():
    t = build_table(w_exact, np.linspace(0.02, 1.0, 99))
    assert t.eta_dec == pytest.approx(0.6, abs=2e-4)          # widths below 1 nm count as no crack
    for e in t.eta[t.eta < t.eta_dec]:
        assert t(float(e)) == pytest.approx(w_exact(float(e)), rel=1e-12)
    assert t(0.61) == 0.0 and t(1.0) == 0.0
    assert t.is_nonincreasing()


def test_interpolation_error_falls_with_the_grid():
    errs = []
    for n in (12, 25, 99):
        t = build_table(w_exact, np.linspace(0.02, 1.0, n))
        knots = t.eta[t.eta < t.eta_dec]
        errs.append(max_rel_error(t, w_exact, 0.5 * (knots[:-1] + knots[1:])))
    assert errs[0] > errs[1] > errs[2] and errs[2] < 5e-3


def test_monte_carlo_on_the_table_matches_the_exact_map():
    t = build_table(w_exact, np.linspace(0.02, 1.0, 99))
    a = monte_carlo_eta(t, 0.15e-3, 0.04e-3, n=2000).summary()
    b = monte_carlo_eta(w_exact, 0.15e-3, 0.04e-3, n=2000).summary()
    for q in ("p05", "p50", "p95"):
        assert a[q] == pytest.approx(b[q], abs=2e-3)


def test_bad_forward_maps_are_refused():
    with pytest.raises(ValueError, match="not finite"):
        build_table(lambda e: float("nan") if e < 0.3 else w_exact(e), np.linspace(0.02, 1.0, 20))
    with pytest.raises(ValueError, match="ill posed"):
        build_table(lambda e: 1e-3 * max(0.0, 0.6 - e) * (1.0 + 5.0 * (e > 0.3)), np.linspace(0.02, 1.0, 20))
    with pytest.raises(ValueError):
        build_table(w_exact, [0.5, 0.2])
