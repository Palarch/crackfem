"""secant_iterate: the generic secant loop (any FE model) that returns non-convergence instead
of raising. Checked against the Stage-1 solve_secant on the same beams."""
import numpy as np
import pytest

from crackfem.model import BeamModel, simply_supported, with_loads
from crackfem.nonlinear import secant_iterate, solve_secant
from crackfem.solve import midpoint_moments, solve

L, EI0, W = 24.5, 7.5e9, 28e3
MREF = W * L**2 / 8


def softening(M):
    return EI0 / (1.0 + np.abs(M) / MREF)


def beam_analyse(model):
    def run(EI):
        sol = solve(model, EI)
        return sol, midpoint_moments(model, sol)
    return run


def test_determinate_beam_one_pass():
    m = with_loads(simply_supported(L, 20), q=-W)
    res = secant_iterate(beam_analyse(m), softening, np.full(m.n_el, EI0))
    assert res.converged and res.reason == "converged"
    assert res.iterations == 1


def test_matches_solve_secant_on_an_indeterminate_beam():
    n = 20
    m = with_loads(BeamModel(x=np.linspace(0, L, n + 1), fixed_dofs=(0, 1, 2 * n)), q=-W)
    ref = solve_secant(m, softening, EI0, tol=1e-10)
    res = secant_iterate(beam_analyse(m), softening, np.full(n, EI0), tol=1e-10)
    assert res.converged and res.iterations == ref.iterations >= 2
    np.testing.assert_allclose(res.EI, ref.EI, rtol=1e-12)
    np.testing.assert_allclose(res.history, ref.history, rtol=1e-9, atol=1e-15)
    assert res.ei_residual[-1] <= 1e-6


def test_two_cycle_is_returned_as_not_converged_with_its_history():
    n = 40
    m = with_loads(BeamModel(x=np.linspace(0, 2 * L, n + 1), fixed_dofs=(0, n, 2 * n)), q=-W)
    Mc = 0.75 * np.max(np.abs(midpoint_moments(m, solve(m, EI0))))
    step = lambda M: np.where(np.abs(M) < Mc, EI0, EI0 / 4)
    res = secant_iterate(beam_analyse(m), step, np.full(n, EI0), max_iter=30)
    assert not res.converged and res.reason == "max_iter"
    assert res.iterations == 30 and len(res.history) == 30 == len(res.ei_residual)


def test_under_relaxation_reaches_the_same_fixed_point():
    n = 20
    m = with_loads(BeamModel(x=np.linspace(0, L, n + 1), fixed_dofs=(0, 1, 2 * n)), q=-W)
    a = secant_iterate(beam_analyse(m), softening, np.full(n, EI0), tol=1e-10)
    b = secant_iterate(beam_analyse(m), softening, np.full(n, EI0), tol=1e-10, relax=0.5)
    assert b.converged and b.iterations > a.iterations
    np.testing.assert_allclose(b.EI, a.EI, rtol=1e-5)          # both stop at tol_EI = 1e-6
    with pytest.raises(ValueError):
        secant_iterate(beam_analyse(m), softening, np.full(n, EI0), relax=0.0)
