"""Grid assembly and partition solve: Logan Example 5.5 through the whole path, and the
reductions a grid must obey (a single line is a beam, a cantilever twists by TL/GJ)."""
import json
from pathlib import Path

import numpy as np
import pytest

from crackfem.gridmodel import GridModel, with_loads
from crackfem.gridsolve import line_node_moments, midpoint_moments, moment_at, solve, stiffness, torque
from crackfem.model import simply_supported
from crackfem.model import with_loads as beam_loads
from crackfem.params import FT, INCH, KIP, KSI
from crackfem.solve import node_moments
from crackfem.solve import solve as beam_solve

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))
L, EI, GJ, W = 24.5, 7.5e9, 1.2e9, 28e3


def test_logan_example_5_5_through_assembly_and_partition():
    ex = FIX["logan_grid_example"]
    names = sorted(ex["nodes_ft"], key=int)                       # "1".."4" -> indices 0..3
    xz = np.array([ex["nodes_ft"][n] for n in names]) * FT
    members = np.array(ex["elements"]) - 1
    fixed = tuple(3 * (n - 1) + d for n in ex["fixed_nodes"] for d in range(3))
    model = with_loads(GridModel(xz, members, fixed), nodal_loads=[(0, ex["load"]["F_y_kip"] * KIP)])
    E, G = ex["E_ksi"] * KSI, ex["G_ksi"] * KSI
    sol = solve(model, E * ex["I_in4"] * INCH**4, G * ex["J_in4"] * INCH**4)
    p = ex["displacements"]
    np.testing.assert_allclose([sol.D[0] / INCH, sol.D[1], sol.D[2]], [p["d1y_in"], p["phi1x_rad"], p["phi1z_rad"]], rtol=5e-3)
    to_us = np.array([KIP, KIP * INCH, KIP * INCH, KIP, KIP * INCH, KIP * INCH])
    for e in range(3):
        np.testing.assert_allclose(sol.end_forces[e] / to_us, ex["local_forces_kip_kin"][str(e + 1)], rtol=5e-3)
    assert sol.R[0::3].sum() == pytest.approx(-ex["load"]["F_y_kip"] * KIP, rel=1e-12)


def _line(n, along_z=False):
    s = np.linspace(0.0, L, n + 1)
    xz = np.stack([np.zeros_like(s), s], axis=1) if along_z else np.stack([s, np.zeros_like(s)], axis=1)
    members = np.array([(k, k + 1) for k in range(n)])
    twist = 2 if along_z else 1                                  # the DOF about the line's own axis
    return GridModel(xz, members, (0, twist, 3 * n))            # w at both ends, twist at one end


@pytest.mark.parametrize("along_z", [False, True])
def test_a_single_line_is_the_beam(along_z):
    n, pts = 8, [(5.0, -80e3), (17.0, -40e3)]
    beam = beam_loads(simply_supported(L, n), q=-W, point_loads=pts)
    ref = beam_solve(beam, EI)
    grid = _line(n, along_z)
    gl = [(int(xp // (L / n)), xp - (L / n) * int(xp // (L / n)), P) for xp, P in pts]
    grid = with_loads(grid, q=-W, point_loads=gl)
    sol = solve(grid, EI, GJ)
    np.testing.assert_allclose(sol.D[0::3], ref.D[0::2], rtol=1e-10, atol=1e-14)
    # bending rotation: theta_z for a line along x; along z the slope dw/dz is -theta_x
    rot = -sol.D[1::3] if along_z else sol.D[2::3]
    np.testing.assert_allclose(rot, ref.D[1::2], rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(line_node_moments(grid, sol, range(n)), node_moments(beam, ref), rtol=1e-9, atol=1e-3)
    np.testing.assert_allclose(midpoint_moments(grid, sol, range(n)),
                               [moment_at(grid, sol, m, L / n / 2) for m in range(n)], rtol=1e-12)


def test_cantilever_twists_by_TL_over_GJ():
    n, T = 4, 15e3
    s = np.linspace(0.0, 3.0, n + 1)
    model = GridModel(np.stack([s, np.zeros_like(s)], axis=1), np.array([(k, k + 1) for k in range(n)]), (0, 1, 2))
    sol = solve(with_loads(model, nodal_loads=[(3 * n + 1, T)]), EI, GJ)
    assert sol.D[3 * n + 1] == pytest.approx(T * 3.0 / GJ, rel=1e-12)
    assert all(torque(sol, m) == pytest.approx(T, rel=1e-12) for m in range(n))
    assert np.allclose(sol.D[0::3], 0.0)                                  # pure twist, no bending


def _two_by_two():
    xz = np.array([(0, 0), (6, 0), (0, 4), (6, 4)], dtype=float)
    members = np.array([(0, 1), (2, 3), (0, 2), (1, 3)])
    fixed = (0, 3, 6)                                                    # three corners vertically
    return with_loads(GridModel(xz, members, fixed + (1, 2)), q=[-5e3, -2e3, 0.0, -1e3],
                      point_loads=[(0, 2.0, -30e3)], nodal_loads=[(9, -10e3), (10, 4e3)])


def _skew_cantilever():
    """Reviewer's probe: a cantilever along (3, 4) m, clamped at the origin, UDL + point load."""
    n = 4
    xz = np.outer(np.linspace(0.0, 1.0, n + 1), [3.0, 4.0])
    return with_loads(GridModel(xz, np.array([(k, k + 1) for k in range(n)]), (0, 1, 2)), q=-2e3,
                      point_loads=[(2, 0.4, -5e3)])


def _resultant(model, sol):
    """Sum of all applied loads and reactions as (F_y, M_x, M_z) about the origin, plus the sum
    of the magnitudes (the scale). A vertical force F_y at (x, z) has the moment r x F =
    (-z F_y, 0, x F_y); member loads are taken as physical loads (UDL resultant at midspan,
    point loads at their position), not their consistent nodal equivalents."""
    C, S = model.cosines
    Lm = model.lengths
    forces = []                                                          # (x, z, F_y)
    couples = []                                                         # (M_x, M_z)
    for m in range(model.n_members):
        i = model.members[m, 0]
        forces.append((*(model.xz[i] + Lm[m] / 2 * np.array([C[m], S[m]])), model.udl()[m] * Lm[m]))
    for m, a, P in model.point_loads:
        forces.append((*(model.xz[model.members[m, 0]] + a * np.array([C[m], S[m]])), P))
    applied = np.zeros(model.n_dof)
    for dof, v in model.nodal_loads:
        applied[dof] += v
    total = applied + sol.R                                              # nodal loads + reactions
    for node in range(model.n_nodes):
        forces.append((*model.xz[node], total[3 * node]))
        couples.append((total[3 * node + 1], total[3 * node + 2]))
    f = np.array(forces)
    c = np.array(couples)
    terms_x = np.concatenate([-f[:, 1] * f[:, 2], c[:, 0]])
    terms_z = np.concatenate([f[:, 0] * f[:, 2], c[:, 1]])
    res = np.array([f[:, 2].sum(), terms_x.sum(), terms_z.sum()])
    scale = np.array([np.abs(f[:, 2]).sum(), np.abs(terms_x).sum(), np.abs(terms_z).sum()])
    return res, scale


def test_two_by_two_grid_is_symmetric_and_in_equilibrium():
    model = _two_by_two()
    st = stiffness(model, [EI, EI, 0.5 * EI, 0.5 * EI], GJ)
    np.testing.assert_allclose(st.K, st.K.T, rtol=1e-12, atol=1e-3)
    sol = solve(model, stiff=st)
    total = 5e3 * 6 + 2e3 * 6 + 1e3 * 4 + 30e3 + 10e3
    assert sol.R[0::3].sum() == pytest.approx(total, rel=1e-12)


@pytest.mark.parametrize("build, EIs", [(_two_by_two, [EI, EI, 0.5 * EI, 0.5 * EI]), (_skew_cantilever, 7e9)])
def test_global_force_and_moment_equilibrium(build, EIs):
    """Reaction forces and moments (at the w, theta_x, theta_z DOFs) balance every applied load,
    moments taken about the origin as r x F."""
    model = build()
    sol = solve(model, EIs, GJ)
    res, scale = _resultant(model, sol)
    assert np.all(np.abs(res) <= 1e-9 * scale), (res, scale)
    assert np.all(scale[1:] > 0.0)

def test_reused_factorisation_gives_the_same_answer():
    m = with_loads(_line(6), q=-W)
    st = stiffness(m, EI, GJ)
    a, b = solve(m, EI, GJ), solve(m, stiff=st)
    np.testing.assert_array_equal(a.D, b.D)


def test_invalid_input_is_refused():
    with pytest.raises(ValueError):
        GridModel(np.zeros((2, 2)), np.array([(0, 1)]), (0,))            # zero-length member
    with pytest.raises(ValueError):
        solve(with_loads(_line(4), q=-W), -1.0, GJ)
    with pytest.raises(ValueError):
        with_loads(_line(4), point_loads=[(0, 99.0, -1.0)])


def test_solve_needs_exactly_one_stiffness_source_of_this_model():
    m = with_loads(_line(4), q=-W)
    st = stiffness(m, EI, GJ)
    with pytest.raises(ValueError, match="give EI and GJ"):
        solve(m)                                                         # neither
    with pytest.raises(ValueError, match="give EI and GJ"):
        solve(m, EI)                                                     # GJ missing
    with pytest.raises(ValueError, match="not both"):
        solve(m, EI, GJ, stiff=st)                                       # both
    with pytest.raises(ValueError, match="does not belong"):
        solve(with_loads(_line(6), q=-W), stiff=st)                      # another model's K


def test_a_mechanism_is_refused_as_not_positive_definite():
    s = np.linspace(0.0, L, 3)
    free_twist = GridModel(np.stack([s, np.zeros_like(s)], axis=1), np.array([(0, 1), (1, 2)]), (0, 6))
    with pytest.raises(ValueError, match="not positive definite"):
        stiffness(free_twist, EI, GJ)                                    # no support resists twist


def test_loads_given_to_the_constructor_obey_the_with_loads_rules():
    base = _line(4)
    direct = GridModel(base.xz, base.members, base.fixed_dofs, q=-W, point_loads=[(1, 2.0, -5e3)],
                       nodal_loads=[(4, 1e3)])
    via = with_loads(base, q=-W, point_loads=[(1, 2.0, -5e3)], nodal_loads=[(4, 1e3)])
    np.testing.assert_array_equal(direct.q, np.full(4, -W))
    assert direct.point_loads == via.point_loads == ((1, 2.0, -5e3),)
    assert direct.nodal_loads == via.nodal_loads == ((4, 1e3),)
    np.testing.assert_array_equal(solve(direct, EI, GJ).D, solve(via, EI, GJ).D)
    with pytest.raises(ValueError):
        GridModel(base.xz, base.members, base.fixed_dofs, point_loads=[(0, 99.0, -1.0)])   # off the member
    with pytest.raises(ValueError):
        GridModel(base.xz, base.members, base.fixed_dofs, point_loads=[(7, 1.0, -1.0)])    # no member 7
    with pytest.raises(ValueError):
        GridModel(base.xz, base.members, base.fixed_dofs, nodal_loads=[(99, 1.0)])         # no DOF 99
    with pytest.raises(ValueError):
        GridModel(base.xz, base.members, base.fixed_dofs, q=[1.0, 2.0])                    # wrong length


def test_line_node_moments_needs_an_ordered_chain():
    m = with_loads(_line(4), q=-W)
    sol = solve(m, EI, GJ)
    with pytest.raises(ValueError, match="same-direction chain"):
        line_node_moments(m, sol, [0, 2, 3])                             # gap
    with pytest.raises(ValueError, match="same-direction chain"):
        line_node_moments(m, sol, [1, 0])                                # reversed order
