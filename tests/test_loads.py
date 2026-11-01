import pytest

from crackfem.loads import (axle_offsets, dead_load_line, impact_factor, max_truck_moment_at,
                            midspan_moment_statics, service_truck_loads, stage_a_line, stage_b_line,
                            truck_fraction, truck_position_max_midspan)
from crackfem.params import FT, KIP, Inputs

INP = Inputs()


def _il(L, x, s):
    """Hand influence line for the sagging moment at x of a simply supported span."""
    return s * (L - x) / L if s <= x else x * (L - s) / L


def test_axle_loads_in_newtons():
    assert [P for _, P in axle_offsets(INP)] == pytest.approx([8 * KIP, 32 * KIP, 32 * KIP])
    assert axle_offsets(INP, reverse=True)[0][1] == pytest.approx(32 * KIP)


def test_axle_offsets_are_0_14_and_28_ft():
    assert [o for o, _ in axle_offsets(INP)] == pytest.approx([0.0, 14 * FT, 28 * FT])


def test_impact_factor_for_s01():
    assert impact_factor(INP) == pytest.approx(50 / (24.5 / FT + 125))


def test_truck_fraction_is_half_of_s_over_5_5():
    assert truck_fraction(INP) == pytest.approx((2.2 / FT) / 5.5 / 2)


def test_scan_finds_the_apex_placement_maximum():
    L = INP.span
    il = lambda x: x / 2 if x <= L / 2 else (L - x) / 2
    best = 0.0
    for rev in (False, True):
        axles = axle_offsets(INP, rev)
        for o_k, _ in axles:
            x0 = L / 2 - o_k
            best = max(best, sum(P * il(x0 + o) for o, P in axles if 0 <= x0 + o <= L))
    _, m = truck_position_max_midspan(INP)
    assert m == pytest.approx(best, rel=1e-6)
    # hand influence line: 32 kip at the apex (6.125 m) + 8 and 32 kip at 14 ft either side (3.9914 m each) = 355.656 kip.m
    assert m == pytest.approx(1_582_036.7, rel=1e-7)


def test_dead_load_is_stage_a_plus_stage_b():
    A = 0.5
    a = INP.unit_weight_concrete * (A + 0.2 * 2.2)
    b = INP.unit_weight_wearing * 0.05 * 2.2 + 2 * INP.barrier_weight / 5
    assert stage_a_line(INP, A) == pytest.approx(a)
    assert stage_b_line(INP) == pytest.approx(b)
    assert dead_load_line(INP, A) == pytest.approx(a + b)


def test_service_moment_combines_dead_and_scaled_truck():
    pos, m_truck = truck_position_max_midspan(INP)
    scale = truck_fraction(INP) * (1 + impact_factor(INP)) * INP.live_load_factor
    M = midspan_moment_statics(INP.span, 20e3, service_truck_loads(INP))
    assert M == pytest.approx(20e3 * INP.span**2 / 8 + scale * m_truck, rel=1e-9)


def test_truck_envelope_at_midspan_equals_the_midspan_scan():
    assert max_truck_moment_at(INP, INP.span / 2) == pytest.approx(truck_position_max_midspan(INP)[1], rel=1e-9)


@pytest.mark.parametrize("f", [0.05, 0.2, 0.3, 0.4, 0.45, 0.6, 0.8, 0.95])
def test_truck_envelope_is_at_least_the_midspan_placement(f):
    L = INP.span
    x = f * L
    pos, _ = truck_position_max_midspan(INP)
    assert max_truck_moment_at(INP, x) >= sum(P * _il(L, x, s) for s, P in pos) * (1 - 1e-12)


@pytest.mark.parametrize("f", [0.1, 0.25, 0.4, 0.45])
def test_truck_envelope_is_symmetric(f):
    L = INP.span
    assert max_truck_moment_at(INP, f * L) == pytest.approx(max_truck_moment_at(INP, (1 - f) * L), rel=1e-9)


@pytest.mark.parametrize("f", [0.1, 0.25, 0.4])
def test_truck_envelope_equals_the_best_placement_with_an_axle_at_the_section(f):
    # moment at x is piecewise linear in truck position; its concave kinks are where an axle crosses x
    L = INP.span
    x = f * L
    best = 0.0
    for rev in (False, True):
        axles = axle_offsets(INP, rev)
        for o_k, _ in axles:
            x0 = x - o_k
            best = max(best, sum(P * _il(L, x, x0 + o) for o, P in axles if 0 <= x0 + o <= L))
    assert max_truck_moment_at(INP, x) == pytest.approx(best, rel=1e-9)
