"""The secant loop on the indeterminate grillage (spec P5): more than one pass, consistent EI,
load shed (not lost), sticky pre-cracked flags, w_pred and the inverse through the grillage."""
import numpy as np
import pytest

from crackfem.grillage import Grillage
from crackfem.inverse import estimate_eta, is_monotone_nonincreasing


@pytest.fixture(scope="module")
def G():
    return Grillage(n_el=10)


@pytest.fixture(scope="module")
def G20():
    return Grillage(n_el=20, n_strands=(30, 29, 27, 29, 30))     # the default design, no design rounds


@pytest.fixture(scope="module")
def damaged(G):
    return G.analyse(G.service_model(0), G.etas(0, 0.3))


def test_loop_needs_more_than_one_pass_and_ends_with_consistent_stiffness(G, damaged):
    r = damaged.result
    assert r.converged and r.iterations >= 2
    assert r.history[-1] <= 1e-6 and r.ei_residual[-1] <= 1e-6
    fn = G.stiffness_fn(G.etas(0, 0.3), damaged.precracked)      # the stiffness belonging to the final moments
    np.testing.assert_allclose(fn(r.M), r.EI, rtol=1e-5)


def test_cracked_girder_sheds_moment_that_the_others_pick_up(G, damaged):
    intact = G.analyse(G.service_model(0), (1.0,) * 5)
    assert damaged.M_B_mid[0] < intact.M_B_mid[0]
    assert damaged.M_B_mid[1] > intact.M_B_mid[1]
    # load is shed, not lost: the girders' midspan moments still add up to about the same
    assert damaged.M_B_mid.sum() == pytest.approx(intact.M_B_mid.sum(), rel=1e-2)
    assert damaged.deflection_mid[0] > intact.deflection_mid[0]


def test_sticky_flag_keeps_an_element_precracked(G):
    gd = G.girders[2]
    etas = G.etas(2, 0.45)
    sound, dec = G._mcr(2, 0.45)
    i = G.deck.mid                                         # element just right of midspan
    n = G.n_el
    assert dec[i] < sound[i]
    M_b = np.zeros(5 * n)
    fn = G.stiffness_fn(etas, np.zeros((5, n), dtype=bool))
    M_b[2 * n + i] = 1.05 * sound[i] - gd.M_A_el[i]        # above the sound Mcr: flag set
    fn(M_b)
    M_b[2 * n + i] = 0.5 * (sound[i] + dec[i]) - gd.M_A_el[i]   # back between decompression and Mcr
    sticky = fn(M_b)[2 * n + i]
    fresh = G.stiffness_fn(etas, np.zeros((5, n), dtype=bool))(M_b)[2 * n + i]
    assert fresh == pytest.approx(gd.EI_I) and sticky < 0.99 * gd.EI_I
    assert fn.precracked[2, i]


def test_w_pred_through_the_grillage(G):
    assert G.w_pred(1.0, 0) == 0.0
    assert G.w_pred(0.2, 0) > 0.0
    assert is_monotone_nonincreasing(G.w_pred_fn(2), np.linspace(0.05, 1.0, 12))
    assert G.w_pred(0.3, 0) == pytest.approx(G.w_pred(0.3, 4), rel=1e-9)    # G1 and G5 are mirror images
    d = G.w_pred_detail(0.3, 0)
    assert d["converged"] and d["iterations"] >= 2


def test_inverse_round_trip_through_the_grillage(G):
    est = estimate_eta(G.w_pred_fn(2), G.w_pred(0.33, 2))
    assert est.status == "ok"
    assert est.eta == pytest.approx(0.33, abs=1e-5)


def test_non_convergence_is_reported_not_raised(G):
    an = G.analyse(G.service_model(0), G.etas(0, 0.3), max_iter=1)
    assert not an.result.converged and an.result.reason == "max_iter" and len(an.result.history) == 1


def test_torsion_sensitivity_and_indicator(G):
    base = G.indicator(0.3, 0)
    soft = G.indicator(0.3, 0, gj_factor=0.5)
    assert base["converged"] == {"dead": True, "service": True} and soft["converged"]["service"]
    assert soft["M_B_service_mid_kNm"]["damaged"][0] != pytest.approx(base["M_B_service_mid_kNm"]["damaged"][0], rel=1e-6)
    assert soft["w_pred_mm"] != pytest.approx(base["w_pred_mm"], rel=1e-6)          # the dead-load solve sees it too
    assert base["cracked_elements"][0] > 0
    assert sum(base["share_service"]["damaged"]) == pytest.approx(1.0)
    assert base["shed_percent_of_Gk"][0] < 0.0 < base["shed_percent_of_Gk"][1]
    assert base["deflection_truck_mm"]["damaged"][0] > base["deflection_truck_mm"]["intact"][0]
    intact = G.indicator(1.0, 0)
    assert intact["cracked_elements"] == [0, 0, 0, 0, 0] and intact["iterations"]["service"] == 1


def test_load_path_is_dead_load_then_equal_truck_steps(G):
    serv = G.service_model(0)
    path = G.load_path(serv, n_steps=4)
    assert len(path) == 5 and path[0] is G.model_dead and path[-1] is serv
    np.testing.assert_allclose(path[2].udl(), G.model_dead.udl())                   # the truck has no UDL
    wheels = lambda m: -sum(P for _, _, P in m.point_loads)                         # total wheel load
    assert wheels(path[2]) == pytest.approx(0.5 * wheels(serv), rel=1e-12)
    assert G.load_path(G.model_dead) == [G.model_dead]                              # at rest: one step
    with pytest.raises(ValueError):
        G.load_path(serv, n_steps=0)


@pytest.mark.parametrize("eta", [0.36, 0.3])
def test_flags_and_moments_do_not_depend_on_the_iteration_path(G, eta):
    a = G.analyse(G.service_model(0), G.etas(0, eta))
    b = G.analyse(G.service_model(0), G.etas(0, eta), relax=0.5)
    assert a.result.converged and b.result.converged
    np.testing.assert_array_equal(a.precracked, b.precracked)
    np.testing.assert_allclose(a.M_B_mid, b.M_B_mid, rtol=1e-5)
    np.testing.assert_allclose(a.result.M, b.result.M, rtol=1e-5, atol=1e-5 * np.abs(a.result.M).max())


@pytest.mark.parametrize("eta", [0.36, 0.3])
def test_every_flag_was_earned_on_the_loading_path(G, eta):
    etas = G.etas(0, eta)
    an = G.analyse(G.service_model(0), etas)
    pre0 = np.array([G.precracked(g, e) for g, e in enumerate(etas)])
    sound = np.array([G._mcr(g, e)[0] for g, e in enumerate(etas)])
    earned = an.path_max > sound
    # a flag not set by the envelope was set by a converged state above the sound Mcr, and back
    np.testing.assert_array_equal(an.precracked & ~pre0, earned & ~pre0)
    assert len(an.step_iterations) == G.n_steps + 1
    # the final state is an equilibrium of the full rule: no unflagged element above its sound Mcr
    total = np.array([gd.M_A_el for gd in G.girders]) + an.result.M.reshape(5, G.n_el)
    assert not np.any((total > sound) & ~an.precracked)


def test_intact_path_is_linear_at_every_step(G):
    an = G.analyse(G.service_model(0), (1.0,) * 5)
    assert an.step_iterations == (1,) * (G.n_steps + 1) and not an.precracked.any()


@pytest.mark.parametrize("eta", [0.36, 0.3])
def test_flags_and_moments_do_not_depend_on_the_load_step_size(G20, eta):
    G = G20                                                       # n_el = 20: where the step size mattered
    a = G.analyse(G.service_model(0), G.etas(0, eta), n_steps=10)
    b = G.analyse(G.service_model(0), G.etas(0, eta), n_steps=40)
    assert a.result.converged and b.result.converged and len(b.step_iterations) == 41
    np.testing.assert_array_equal(a.precracked, b.precracked)
    np.testing.assert_allclose(a.M_B_mid, b.M_B_mid, rtol=1e-5)
    np.testing.assert_allclose(a.result.M, b.result.M, rtol=1e-5, atol=1e-5 * np.abs(a.result.M).max())
