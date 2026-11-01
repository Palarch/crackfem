"""Deck layout, Hambly properties and load mapping: counts, closed-form values and statics."""
import numpy as np
import pytest

from crackfem.deck import (area_to_girders, fixed_properties, girder_torsion_constant, layout,
                           lever_weights, line_loads_model, line_to_girders, n_lanes,
                           stage_b_dead_lines, torsion_constant_open, truck_placements, wheel_loads)
from crackfem.gridsolve import load_vector
from crackfem.materials import ec2_ecm
from crackfem.params import FT, INCH, Inputs

INP = Inputs()


@pytest.fixture(scope="module")
def deck():
    return layout(INP, 10)


def test_layout_counts_and_supports(deck):
    ng, n = INP.n_girders, 10
    assert deck.model.n_nodes == ng * (n + 1)
    assert deck.model.n_members == ng * n + (n + 1) * (ng - 1)
    assert deck.kind.count("diaphragm") == 3 * (ng - 1)
    assert len(deck.model.fixed_dofs) == 2 * ng and all(d % 3 == 0 for d in deck.model.fixed_dofs)
    assert deck.strip_width[deck.cross_members].sum() == pytest.approx(INP.span * (ng - 1))
    np.testing.assert_allclose(deck.girder_widths(), [1.1 + 3 * FT, 2.2, 2.2, 2.2, 1.1 + 3 * FT])
    with pytest.raises(ValueError):
        layout(INP, 7)


def test_type_iv_torsion_constant_is_sum_b_t3_over_3():
    assert girder_torsion_constant(INP) == pytest.approx((26 + 20 + 38) * 8**3 / 3 * INCH**4, rel=1e-12)
    assert torsion_constant_open([(0.2, 1.0)]) == pytest.approx(1.0 * 0.2**3 / 3)


def test_hambly_slab_strip_and_girder_torsion(deck):
    EI, GJ = fixed_properties(deck)
    Es = ec2_ecm(INP.fck_slab)
    Gs = Es / 2.4
    t, b = INP.slab_thickness, INP.span / 10
    m = int(deck.cross_members[3, 1])                       # an interior slab strip
    assert deck.kind[m] == "slab"
    assert EI[m] == pytest.approx(Es * b * t**3 / 12) and GJ[m] == pytest.approx(Gs * b * t**3 / 6)
    d = int(deck.cross_members[5, 1])                       # the midspan diaphragm
    assert deck.kind[d] == "diaphragm"
    assert EI[d] == pytest.approx(Es * (b * t**3 + 0.25 * 1.10**3) / 12)
    assert np.all(np.isnan(EI[deck.girder_members.ravel()]))
    Gg = ec2_ecm(INP.fck_girder) / 2.4
    assert GJ[deck.girder_members[2, 0]] == pytest.approx(Gg * girder_torsion_constant(INP) + Gs * 2.2 * t**3 / 6)


def test_lever_rule_and_overhang_torque():
    z = np.array([0.0, 2.2, 4.4])
    assert lever_weights(z, 0.55) == [(0, 0.75, 0.0), (1, 0.25, 0.0)]
    assert lever_weights(z, -0.5) == [(0, 1.0, -0.5)]
    q, t = line_to_girders(z, -0.5, 10.0)
    assert q.tolist() == [10.0, 0.0, 0.0] and t.tolist() == [-5.0, 0.0, 0.0]


@pytest.mark.parametrize("za, zb", [(-0.6, 9.4), (0.3, 1.9), (-0.9, 0.0), (1.0, 7.0)])
def test_area_load_mapping_keeps_force_and_moment(za, zb):
    z = INP.girder_spacing * np.arange(5)
    w = 1.1e3
    q, t = area_to_girders(z, za, zb, w)
    assert q.sum() == pytest.approx(w * (zb - za), rel=1e-12)
    # moment about the x axis through z = 0: m_x = z P for a downward load P
    assert (q * z).sum() + t.sum() == pytest.approx(w * (zb**2 - za**2) / 2, rel=1e-12, abs=1e-9)


def test_stage_b_dead_load_totals(deck):
    q, t = stage_b_dead_lines(deck)
    lo, hi = deck.curbs
    total = INP.unit_weight_wearing * INP.wearing_thickness * (hi - lo) + INP.n_barriers * INP.barrier_weight
    assert q.sum() == pytest.approx(total, rel=1e-12)
    np.testing.assert_allclose(q, q[::-1], rtol=1e-12)                    # symmetric deck
    assert t[0] < 0.0 and t[-1] == pytest.approx(-t[0]) and np.all(t[1:-1] == 0.0)
    F, _ = load_vector(line_loads_model(deck, q, t))
    assert -F[0::3].sum() == pytest.approx(total * INP.span, rel=1e-12)
    assert F[1::3].sum() == pytest.approx(t.sum() * INP.span, abs=1e-6)


def test_wheel_loads_keep_force_and_moment(deck):
    axles = [(3.0, 100e3), (12.25, 50e3)]
    lines = (-0.3, 1.2, 5.0)                                              # one wheel line on the overhang
    pts, nodal = wheel_loads(deck, axles, lines, 1.5)
    total = 1.5 * sum(P for _, P in axles) / 2 * len(lines)
    assert -sum(P for _, _, P in pts) == pytest.approx(total, rel=1e-12)
    z_of = {int(m): deck.z[g] for g in range(deck.n_girders) for m in deck.girder_members[g]}
    moment = sum(-P * z_of[m] for m, _, P in pts) + sum(v for dof, v in nodal if dof % 3 == 1)
    assert moment == pytest.approx(1.5 * sum(P for _, P in axles) / 2 * sum(lines), rel=1e-12)


def test_truck_placements_respect_curbs_lanes_and_clearance(deck):
    lo, hi = deck.curbs
    assert n_lanes(deck) == 2                                              # 10.0 m roadway
    ps = truck_placements(deck, step=0.1)
    edge, sw, cw = 2 * FT, 6 * FT, 10 * FT
    assert {len(p) for p in ps} == {2, 4}
    for p in ps:
        assert min(p) >= lo + edge - 1e-9 and max(p) <= hi - edge + 1e-9
        for a, b in zip(p[0::2], p[1::2]):
            assert b - a == pytest.approx(sw)
        if len(p) == 4:
            assert p[2] - p[1] >= 2 * edge - 1e-9                          # clearance widths do not overlap
    has = lambda target: any(len(p) == len(target) and np.allclose(p, target, atol=1e-9) for p in ps)
    assert has((lo + edge, lo + edge + sw)) and has((hi - edge - sw, hi - edge))   # both curb positions
    mirrored = {tuple(round(lo + hi - v, 6) for v in reversed(p)) for p in ps}
    assert mirrored == {tuple(round(v, 6) for v in p) for p in ps}
    assert cw == pytest.approx(2 * edge + sw)


def _starts(p, edge):
    return p[0] - edge, p[2] - edge                                         # clearance-width starts


def test_lane_rule_puts_each_truck_in_its_own_12ft_lane(deck):
    lo, hi = deck.curbs
    edge, cw, lane = 2 * FT, 10 * FT, 12 * FT
    lanes = truck_placements(deck, step=0.1)                                # default: lane_rule="lanes"
    assert lanes == truck_placements(deck, step=0.1, lane_rule="lanes")
    two = [p for p in lanes if len(p) == 4]
    assert two
    for p in two:
        s1, s2 = _starts(p, edge)
        # closed form of "two non-overlapping 12 ft lanes between the curbs, each holding its 10 ft truck"
        assert s1 >= lo - 1e-9 and s2 + cw <= hi + 1e-9 and s2 >= s1 + cw - 1e-9
        assert s2 >= lo + lane - 1e-9 and s1 + cw <= hi - lane + 1e-9
    clearance = truck_placements(deck, step=0.1, lane_rule="clearance")
    key = lambda ps: {tuple(round(v, 6) for v in p) for p in ps}
    # off the lane-edge positions (added under "lanes" only), the lane placements are clearance placements
    lane_edge = {round(lo + lane + edge, 6), round(hi - lane - cw + edge, 6)}
    on_grid = {p for p in key(lanes) if len(p) == 2 or not {p[0], p[2]} & lane_edge}
    assert on_grid < key(clearance) and len(key(lanes)) < len(key(clearance))
    assert {p for p in key(clearance) if len(p) == 2} == {p for p in key(lanes) if len(p) == 2}
    # clearance allows a second truck closer to the left curb than one lane width; lanes do not
    assert any(len(p) == 4 and _starts(p, edge)[1] < lo + lane - 1e-6 for p in clearance)
    first_in_lane2 = (lo + edge, lo + edge + 6 * FT, lo + lane + edge, lo + lane + edge + 6 * FT)
    assert any(len(p) == 4 and np.allclose(p, first_in_lane2, atol=1e-9) for p in truck_placements(deck, step=0.05))
    with pytest.raises(ValueError):
        truck_placements(deck, lane_rule="anywhere")


def test_clearance_rule_is_the_previous_placement_set(deck):
    assert len(truck_placements(deck, 0.05, lane_rule="clearance")) == 6757


def test_diaphragm_options_none_and_composite(deck):
    Es = ec2_ecm(INP.fck_slab)
    t, hd, bd = INP.slab_thickness, 1.10, 0.25
    web_EI, web_GJ = fixed_properties(deck)
    assert np.array_equal(web_EI, fixed_properties(deck, diaphragm="web")[0], equal_nan=True)
    mid = [int(m) for m in deck.cross_members[deck.mid]]
    ends = [int(m) for m in deck.cross_members[[0, -1]].ravel()]
    slab = [int(m) for m in deck.cross_members[3]]
    # none: the midspan diaphragm becomes a plain slab strip; the end diaphragms stay
    EI, GJ = fixed_properties(deck, diaphragm="none")
    np.testing.assert_allclose(EI[mid], web_EI[slab], rtol=1e-12)
    np.testing.assert_allclose(GJ[mid], web_GJ[slab], rtol=1e-12)
    np.testing.assert_allclose(EI[ends], web_EI[ends], rtol=1e-12)
    # composite: a T section, flange = the slab strip of its station, web below it
    EI, GJ = fixed_properties(deck, diaphragm="composite")
    for m in mid + ends:
        b = deck.strip_width[m]
        A1, y1, A2, y2 = b * t, hd + t / 2, bd * hd, hd / 2
        yb = (A1 * y1 + A2 * y2) / (A1 + A2)
        I_T = b * t**3 / 12 + A1 * (y1 - yb) ** 2 + bd * hd**3 / 12 + A2 * (y2 - yb) ** 2
        assert EI[m] == pytest.approx(Es * I_T, rel=1e-12) and EI[m] > web_EI[m]
    np.testing.assert_allclose(GJ, web_GJ, rtol=1e-12)
    others = [m for m in deck.cross_members.ravel() if deck.kind[m] == "slab"]
    np.testing.assert_allclose(EI[others], web_EI[others], rtol=1e-12)
    with pytest.raises(ValueError):
        fixed_properties(deck, diaphragm="steel")
