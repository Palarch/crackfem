"""Grillage of S01 at eta = 1: strand design, truck placement, sanity bands, linearity."""
import numpy as np
import pytest

from crackfem.deck import lever_weights
from crackfem.grillage import Grillage
from crackfem.gridsolve import line_node_moments, midpoint_moments, solve
from crackfem.loads import truck_fraction, truck_position_max_midspan
from crackfem.params import FT
from crackfem.staging import prestress_staged


@pytest.fixture(scope="module")
def G():
    return Grillage(n_el=10)


def test_every_girder_has_no_soffit_tension_at_service_and_the_count_is_minimal(G):
    mid = G.deck.mid
    assert all(s <= 0.0 for s in G.soffit_service)
    for g, gd in enumerate(G.girders):
        M_A = float(gd.M_A_nd[mid])
        fewer = prestress_staged(gd.make_section, gd.is_girder, gd.h_girder, G.n_strands[g] - 1, G.inp.strand_area,
                                 G.fpe, M_A, M_A + G.M_B_service_mid_I[g])
        assert fewer.soffit_stress_service > 0.0
    assert G.design_rounds[-1] == G.n_strands


def test_symmetric_deck_gives_symmetric_design_and_moments(G):
    assert G.n_strands == G.n_strands[::-1]
    np.testing.assert_allclose(G.M_B_service_mid_I, G.M_B_service_mid_I[::-1], rtol=1e-9)
    np.testing.assert_allclose(G.envelope_el, G.envelope_el[::-1], rtol=1e-9)


def test_superposed_placement_moments_equal_a_direct_solve(G):
    for p in (0, len(G.placements) // 3, len(G.placements) - 1):
        model = G.truck_model(G.placements[p])
        sol = solve(model, stiff=G.stiff_I)
        direct = [line_node_moments(model, sol, G.deck.girder_members[g])[G.deck.mid] for g in range(G.deck.n_girders)]
        np.testing.assert_allclose(G.M_truck_mid[p], direct, rtol=1e-9)


def test_exterior_girder_is_loaded_from_its_curb(G):
    lo, _ = G.deck.curbs
    lines = G.placements[G.critical[0]]
    assert min(lines) == pytest.approx(lo + 2 * FT)
    assert len(lines) == 4                                                  # two trucks govern G1


def test_truck_share_sanity_bands(G):
    one_truck = truck_position_max_midspan(G.inp)[1] * G.truck_scale
    frac = G.M_truck_mid[[G.critical[g] for g in range(5)], range(5)] / one_truck
    # interior girders: below the AASHTO S/5.5 fraction (a code bound), and not far below it
    assert 0.5 < max(frac[1:4]) / truck_fraction(G.inp) < 1.1
    # exterior girder: between the lever rule (hinged slab) and Courbon (rigid deck), same wheels
    z, lines = G.deck.z, G.placements[G.critical[0]]
    lever = sum(s / 2 for zw in lines for g, s, _ in lever_weights(z, zw) if g == 0)
    zbar = z - z.mean()
    courbon = sum((1 / z.size + (zw - z.mean()) * zbar[0] / (zbar**2).sum()) / 2 for zw in lines)
    assert lever < frac[0] < courbon


def test_intact_analysis_is_the_linear_solve(G):
    model = G.service_model(0)
    an = G.analyse(model, (1.0,) * 5)
    lin = solve(model, G.EI_I, G.GJ0)
    assert an.result.converged and an.result.iterations == 1
    np.testing.assert_allclose(an.result.solution.D, lin.D, rtol=1e-9, atol=1e-12)
    assert not an.precracked.any()


def test_envelope_covers_the_service_moments_and_nothing_is_precracked_at_eta_1(G):
    for g in range(5):
        model = G.truck_model(G.placements[G.critical[g]], base=G.model_dead)
        sol = solve(model, stiff=G.stiff_I)
        service = G.girders[g].M_A_el + midpoint_moments(model, sol, G.deck.girder_members[g])
        assert np.all(G.envelope_el[g] >= service - 1e-6 * np.abs(service).max())
        assert not G.precracked(g, 1.0).any()


def test_default_placements_follow_the_lane_rule_and_clearance_is_more_severe_for_g1(G):
    lo, hi = G.deck.curbs
    lane, cw, edge = 12 * FT, 10 * FT, 2 * FT
    assert G.lane_rule == "lanes"
    for p in G.placements:
        if len(p) == 4:
            assert p[2] - edge >= lo + lane - 1e-9 and p[0] - edge + cw <= hi - lane + 1e-9
    C = Grillage(n_el=10, n_strands=G.n_strands, lane_rule="clearance")
    assert C.M_truck_mid[:, 0].max() > G.M_truck_mid[:, 0].max()


def test_composite_diaphragms_reproduce_the_independent_review_moments():
    # review (independent script, n_el = 20, clearance rule, composite T diaphragms): G1 1838.8, G2 1432.7 kN m
    C = Grillage(n_el=20, n_strands=(31, 29, 27, 29, 31), lane_rule="clearance", diaphragm="composite")
    crit = [C.M_truck_mid[C.critical[g], g] / 1e3 for g in (0, 1)]
    assert crit[0] == pytest.approx(1838.8, rel=5e-3) and crit[1] == pytest.approx(1432.7, rel=5e-3)
    N = Grillage(n_el=10, n_strands=(30, 29, 27, 29, 30), diaphragm="none")
    mid = [int(m) for m in N.deck.cross_members[N.deck.mid]]
    strip = int(N.deck.cross_members[3, 0])
    np.testing.assert_allclose(N.EI_fixed[mid], N.EI_fixed[strip], rtol=1e-12)
