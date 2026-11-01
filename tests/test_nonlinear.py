import numpy as np
import pytest

from crackfem.model import BeamModel, simply_supported, with_loads
from crackfem.nonlinear import solve_secant
from crackfem.solve import midpoint_moments, moment_at, solve

L, EI0, W = 24.5, 7.5e9, 28e3
MREF = W * L**2 / 8


def softening(M):
    return EI0 / (1.0 + np.abs(M) / MREF)


def test_determinate_beam_converges_in_one_pass():
    m = with_loads(simply_supported(L, 20), q=-W)
    res = solve_secant(m, softening, EI0)
    assert res.iterations == 1
    assert res.history[0] < 1e-9   # round-off only; the one-pass property is pinned by iterations == 1


def test_constant_stiffness_reproduces_the_linear_solution():
    m = with_loads(simply_supported(L, 20), q=-W)
    res = solve_secant(m, lambda M: np.full(M.size, EI0), EI0)
    np.testing.assert_allclose(res.solution.D, solve(m, EI0).D, rtol=1e-12)


def test_indeterminate_beam_redistributes_moment_and_reaches_a_fixed_point():
    n = 20
    m = with_loads(BeamModel(x=np.linspace(0, L, n + 1), fixed_dofs=(0, 1, 2 * n)), q=-W)
    res = solve_secant(m, softening, EI0, tol=1e-10)
    assert res.iterations >= 2
    np.testing.assert_allclose(softening(res.M_mid), res.EI, rtol=1e-8)
    assert res.solution.R[0::2].sum() == pytest.approx(W * L, rel=1e-10)   # vertical reactions only (R[1] is the fixed-end moment)
    linear = solve(m, EI0)
    linear_fixed_end = moment_at(m, linear, 0, 0.0)
    ratio = moment_at(m, res.solution, 0, 0.0) / linear_fixed_end
    assert ratio == pytest.approx(0.9663, abs=2e-3)       # the stiffer-in-the-span beam sheds fixed-end moment
    assert midpoint_moments(m, res.solution).max() > midpoint_moments(m, linear).max()   # ...to the sagging span


def test_two_cycle_raises_instead_of_returning_a_non_converged_answer():
    n = 40      # two spans of L, supports at 0, L, 2L
    m = with_loads(BeamModel(x=np.linspace(0, 2 * L, n + 1), fixed_dofs=(0, n, 2 * n)), q=-W)
    Mc = 0.75 * np.max(np.abs(midpoint_moments(m, solve(m, EI0))))      # plain iteration 2-cycles here
    step = lambda M: np.where(np.abs(M) < Mc, EI0, EI0 / 4)
    with pytest.raises(RuntimeError, match="did not converge"):
        solve_secant(m, step, EI0)


def test_stiffness_residual_check_keeps_iterating_after_the_moments_settle():
    # Determinate beam: moments never change, so a moment-only test would stop at once. This stateful
    # stiffness_fn relaxes halfway toward the softening stiffness on every call.
    m = with_loads(simply_supported(L, 20), q=-W)
    state = {"EI": np.full(m.n_el, EI0)}

    def relaxing(M):
        state["EI"] = 0.5 * (state["EI"] + softening(M))
        return state["EI"].copy()

    res = solve_secant(m, relaxing, EI0)
    assert res.iterations > 5
    np.testing.assert_allclose(res.EI, softening(res.M_mid), rtol=1e-5)
