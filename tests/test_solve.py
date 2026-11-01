import numpy as np
import pytest
from scipy.integrate import quad

from crackfem.model import BeamModel, simply_supported, with_loads
from crackfem.solve import midpoint_moments, moment_at, node_moments, solve

L, EI, W = 24.5, 7.5e9, 28e3


@pytest.mark.parametrize("n", [1, 2, 3, 5, 8, 40])
def test_simply_supported_udl_is_exact_at_the_nodes_for_any_mesh(n):
    m = with_loads(simply_supported(L, n), q=-W)
    sol = solve(m, EI)
    x = m.x
    exact = -W * x * (L**3 - 2 * L * x**2 + x**3) / (24 * EI)
    np.testing.assert_allclose(sol.D[0::2], exact, rtol=1e-10, atol=1e-15)


def test_reactions_balance_the_applied_loads():
    m = with_loads(simply_supported(L, 10), q=-W, point_loads=[(3.3, -100e3), (12.25, -50e3)])
    sol = solve(m, EI)
    assert sol.R.sum() == pytest.approx(W * L + 150e3, rel=1e-12)


def test_cantilever_tip_load():
    m = BeamModel(x=np.linspace(0, 3.0, 7), fixed_dofs=(0, 1))
    m = with_loads(m, q=0.0, point_loads=[(3.0, 10e3)])
    sol = solve(m, 2.0e7)
    assert sol.D[-2] == pytest.approx(10e3 * 3.0**3 / (3 * 2.0e7), rel=1e-12)


def test_propped_cantilever_fixed_end_moment():
    m = BeamModel(x=np.linspace(0, L, 9), fixed_dofs=(0, 1, 16))
    m = with_loads(m, q=-W)
    sol = solve(m, EI)
    assert moment_at(m, sol, 0, 0.0) == pytest.approx(-W * L**2 / 8, rel=1e-10)


def test_stepped_stiffness_matches_virtual_work():
    n = 10
    m = with_loads(simply_supported(L, n), q=-W)
    EI_el = np.where(m.midpoints < L / 2, EI, 0.6 * EI)
    sol = solve(m, EI_el)
    M = lambda x: W * x * (L - x) / 2
    unit = lambda x: x / 2 if x <= L / 2 else (L - x) / 2
    d1 = quad(lambda x: M(x) * unit(x) / EI, 0, L / 2)[0]
    d2 = quad(lambda x: M(x) * unit(x) / (0.6 * EI), L / 2, L)[0]
    assert -sol.D[2 * (n // 2)] == pytest.approx(d1 + d2, rel=1e-9)


def test_midpoint_moments_match_statics_with_point_loads():
    pts = [(5.0, -80e3), (17.0, -40e3)]
    m = with_loads(simply_supported(L, 12), q=-W, point_loads=pts)
    sol = solve(m, EI)
    RA = W * L / 2 + sum(-P * (L - xp) / L for xp, P in pts)
    def statics(x):
        return RA * x - W * x * x / 2 - sum(-P * (x - xp) for xp, P in pts if xp < x)
    np.testing.assert_allclose(midpoint_moments(m, sol), [statics(x) for x in m.midpoints], rtol=1e-9)


def test_node_moments_match_statics():
    pts = [(5.0, -80e3), (17.0, -40e3)]
    m = with_loads(simply_supported(L, 12), q=-W, point_loads=pts)
    sol = solve(m, EI)
    RA = W * L / 2 + sum(-P * (L - xp) / L for xp, P in pts)
    statics = [RA * x - W * x * x / 2 - sum(-P * (x - xp) for xp, P in pts if xp < x) for x in m.x]
    np.testing.assert_allclose(node_moments(m, sol), statics, rtol=1e-9, atol=1e-6)


def _statics(xs, pts):
    RA = W * L / 2 + sum(-P * (L - xp) / L for xp, P in pts)
    return np.array([RA * x - W * x * x / 2 - sum(-P * (x - xp) for xp, P in pts if xp < x) for x in xs])


def test_point_loads_on_nodes_and_ends_with_varying_EI():
    n = 12
    m0 = simply_supported(L, n)
    pts = [(0.0, -10e3), (m0.x[4], -80e3), (m0.x[7], -40e3), (L, -10e3)]
    m = with_loads(m0, q=-W, point_loads=pts)
    sol = solve(m, EI * np.linspace(1.0, 0.4, n))          # determinate: M must not depend on EI
    np.testing.assert_allclose(node_moments(m, sol), _statics(m.x, pts), rtol=1e-9, atol=1e-5)
    np.testing.assert_allclose(midpoint_moments(m, sol), _statics(m.midpoints, pts), rtol=1e-9)


def test_moment_is_continuous_and_nodal_equilibrium_holds():
    m = with_loads(simply_supported(L, 12), q=-W, point_loads=[(5.0, -80e3)])
    ef = solve(m, EI).end_forces
    np.testing.assert_allclose(ef[:-1, 3] + ef[1:, 1], 0.0, atol=1e-5)   # M2(e-1) = -M1(e)


def test_fixed_support_reactions_force_and_moment():
    n = 6
    m = with_loads(BeamModel(x=np.linspace(0, 3.0, n + 1), fixed_dofs=(0, 1)), q=-5e3, point_loads=[(3.0, -10e3)])
    sol = solve(m, 2e7)
    assert sol.R[0] == pytest.approx(5e3 * 3.0 + 10e3, rel=1e-12)
    assert sol.R[1] == pytest.approx(5e3 * 4.5 + 10e3 * 3.0, rel=1e-12)  # CCW on the beam


def test_non_increasing_node_coordinates_are_rejected():
    for x in ([0.0, 1.0, 1.0, 2.0], [0.0, 2.0, 1.0], [0.0]):
        with pytest.raises(ValueError, match="strictly increasing"):
            BeamModel(x=np.array(x), fixed_dofs=(0, 2))


@pytest.mark.parametrize("fixed", [(0,), (), (0, 99), (-1, 2)])
def test_fewer_than_two_or_out_of_range_restraints_are_rejected(fixed):
    with pytest.raises(ValueError, match="fixed_dofs"):
        BeamModel(x=np.linspace(0.0, 3.0, 4), fixed_dofs=fixed)


@pytest.mark.parametrize("bad", [-1.0, 0.0, np.nan, np.inf])
def test_non_positive_or_non_finite_EI_is_rejected(bad):
    m = with_loads(simply_supported(3.0, 4), q=-1e3)
    with pytest.raises(ValueError, match="EI must be finite and positive"):
        solve(m, bad)
    with pytest.raises(ValueError, match="EI must be finite and positive"):
        solve(m, np.array([1e7, 1e7, bad, 1e7]))


def test_point_load_within_rounding_of_an_end_is_clipped_but_beyond_is_rejected():
    base = simply_supported(L, 8)
    at_end = solve(with_loads(base, point_loads=[(L, -50e3)]), EI)
    just_past = solve(with_loads(base, point_loads=[(L + 1e-12, -50e3)]), EI)
    np.testing.assert_array_equal(just_past.D, at_end.D)
    np.testing.assert_array_equal(just_past.end_forces, at_end.end_forces)
    at_start = solve(with_loads(base, point_loads=[(0.0, -50e3)]), EI)
    just_before = solve(with_loads(base, point_loads=[(-1e-12, -50e3)]), EI)
    np.testing.assert_array_equal(just_before.D, at_start.D)
    for bad in (L + 1e-6, -1e-6):
        with pytest.raises(ValueError, match="outside the beam"):
            solve(with_loads(base, point_loads=[(bad, -50e3)]), EI)
