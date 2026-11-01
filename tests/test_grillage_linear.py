"""Linear grillage checks against independent closed forms (spec §6, grillage rows):
identical girders share equally, a lone girder line is the Stage-1 girder, and a torsion-free
deck between flexible and rigid cross members reproduces the lever rule and Courbon's rule."""
import numpy as np
import pytest

from crackfem.deck import fixed_properties, layout
from crackfem.girder1d import Girder1D
from crackfem.gridmodel import GridModel, with_loads
from crackfem.gridsolve import line_node_moments, midpoint_moments, solve, stiffness
from crackfem.nonlinear import secant_iterate
from crackfem.params import Inputs

INP = Inputs()
N = 10
EI_G = 2.6e10          # any girder EI; the checks are ratios or closed forms in EI


@pytest.fixture(scope="module")
def deck():
    return layout(INP, N)


def _props(deck, cross_scale=1.0, gj_scale=1.0):
    EI, GJ = fixed_properties(deck)
    EI = EI * cross_scale
    EI[deck.girder_members.ravel()] = EI_G
    return EI, GJ * gj_scale


def test_identical_girders_under_equal_load_share_equally(deck):
    q = 20e3
    gm = deck.girder_members.ravel()
    qm = np.zeros(deck.model.n_members)
    qm[gm] = -q
    model = with_loads(deck.model, q=qm)
    EI, GJ = _props(deck)
    st = stiffness(model, EI, GJ)
    np.testing.assert_allclose(st.K, st.K.T, rtol=1e-10, atol=1e-2)
    sol = solve(model, stiff=st)
    L = INP.span
    for g in range(deck.n_girders):
        m = line_node_moments(model, sol, deck.girder_members[g])[deck.mid]
        assert m == pytest.approx(q * L**2 / 8, rel=1e-9)
        assert -sol.D[3 * deck.node(g, deck.mid)] == pytest.approx(5 * q * L**4 / (384 * EI_G), rel=1e-9)
    assert sol.R[0::3].sum() == pytest.approx(deck.n_girders * q * L, rel=1e-12)


def test_a_lone_girder_line_is_the_stage1_girder():
    g = Girder1D(n_el=N)
    x = g.model_service.x
    line = GridModel(np.stack([x, np.zeros_like(x)], axis=1), np.array([(k, k + 1) for k in range(N)]),
                     (0, 1, 3 * N))
    dx = x[1] - x[0]
    pts = [(min(int(xp // dx), N - 1), xp - dx * min(int(xp // dx), N - 1), P) for xp, P in g.model_service.point_loads]
    model = with_loads(line, q=-g.q_b, point_loads=pts)
    GJ = 1.0e9
    for eta in (1.0, 0.36):
        ref = g.analyse(eta, g.model_service)

        def run(EI):
            sol = solve(model, EI, GJ)
            return sol, midpoint_moments(model, sol, range(N))

        res = secant_iterate(run, g.stiffness_fn(eta, g.precracked_zone(eta)), np.full(N, g.EI_I))
        assert res.converged and res.iterations == ref.iterations == 1
        np.testing.assert_allclose(res.EI, ref.EI, rtol=1e-9)
        np.testing.assert_allclose(res.solution.D[0::3], ref.solution.D[0::2], rtol=1e-9, atol=1e-12)


def _girder_shares(deck, model, EI, GJ):
    sol = solve(model, EI, GJ)
    m = np.array([line_node_moments(model, sol, deck.girder_members[g])[deck.mid] for g in range(deck.n_girders)])
    return m / m.sum()


def test_torsion_free_deck_between_lever_rule_and_courbon(deck):
    P = 100e3
    model = with_loads(deck.model, nodal_loads=[(3 * deck.node(0, deck.mid), -P)])   # on G1 at midspan
    zbar = deck.z - deck.z.mean()
    courbon = 1 / deck.n_girders + zbar[0] * zbar / (zbar**2).sum()
    rigid = _girder_shares(deck, model, *_props(deck, cross_scale=1e4, gj_scale=1e-6))
    np.testing.assert_allclose(rigid, courbon, atol=5e-3)
    flexible = _girder_shares(deck, model, *_props(deck, cross_scale=1e-8, gj_scale=1e-6))
    assert flexible[0] == pytest.approx(1.0, abs=1e-3)                      # lever rule: all on G1
    real = _girder_shares(deck, model, *_props(deck))
    assert courbon[0] < real[0] < flexible[0]                               # a real deck lies between
