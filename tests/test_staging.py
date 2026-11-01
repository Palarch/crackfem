import numpy as np
import pytest

from crackfem.geometry import discretize, gross_properties, slab
from crackfem.section import (PCSection, concrete_stress, resultants, soffit_stress, solve_state,
                              uncracked_EI)
from crackfem.staging import design_prestress_staged, girder_only, staged_section

E_G, E_S, E_P = 30e9, 27e9, 195e9
H_G, Y_P, A_P = 0.8, 0.1, 1.5e-3
GIRDER = slab(0.0, H_G, 0.4, E=E_G, fct=3e6)        # a rectangle is a trapezoid with b0 == b1
DECK = slab(H_G, 0.15, 1.2, E=E_S, fct=2.5e6)
FIB = discretize([GIRDER, DECK], 2e-4)
IS_GIRDER = FIB.y < H_G
H = H_G + 0.15
E0, M_A, DM = 4e-3, 200e3, 300e3


def composite(A_p=A_P):
    return PCSection(FIB, Y_P, A_p, E_P, H)


def transformed(parts):
    """Transformed section in girder concrete (n = E / E_G), strand included."""
    A, y, I = gross_properties(parts, E_ref=E_G)
    n = E_P / E_G
    At = A + n * A_P
    yt = (A * y + n * A_P * Y_P) / At
    It = I + A * (y - yt) ** 2 + n * A_P * (Y_P - yt) ** 2
    return At, yt, It


def test_composite_solve_at_the_stage_a_moment_reproduces_stage_a_with_an_unstressed_slab():
    sec = staged_section(composite(), IS_GIRDER, H_G, E0, M_A)
    a = solve_state(girder_only(composite(), IS_GIRDER, H_G), E0, M_A, cracked=False)
    b = solve_state(sec, E0, M_A, cracked=False)
    assert (b.eps0, b.kappa) == pytest.approx((a.eps0, a.kappa), rel=1e-9)
    assert np.abs(concrete_stress(sec, b)[~IS_GIRDER]).max() < 1.0      # Pa


def test_stage_b_increment_matches_the_two_stage_hand_calculation():
    At_g, yt_g, It_g = transformed([GIRDER])
    At_c, yt_c, It_c = transformed([GIRDER, DECK])
    P_lock = E_P * A_P * E0                     # force the locked-in strand strain applies (not Ch3's {f0})
    e_g = yt_g - Y_P
    soffit_a = -P_lock / At_g - P_lock * e_g * yt_g / It_g + M_A * yt_g / It_g     # stage A, bare girder
    soffit_b = DM * yt_c / It_c                                          # stage B, composite
    sec = staged_section(composite(), IS_GIRDER, H_G, E0, M_A)
    a = solve_state(sec, E0, M_A, cracked=False)
    b = solve_state(sec, E0, M_A + DM, cracked=False)
    assert soffit_stress(sec, b) == pytest.approx(soffit_a + soffit_b, rel=1e-5)
    assert b.kappa - a.kappa == pytest.approx(DM / (E_G * It_c), rel=1e-5)
    assert uncracked_EI(sec) == pytest.approx(E_G * It_c, rel=1e-5)


@pytest.mark.parametrize("M", [1.2e6, 1.8e6])
def test_cracked_staged_state_is_in_equilibrium_with_no_concrete_tension(M):
    sec = staged_section(composite(), IS_GIRDER, H_G, E0, M_A)
    st = solve_state(sec, E0, M, cracked=True)
    N, M_int = resultants(sec, E0, st)
    assert abs(N) < 1e-6 * E_P * A_P * E0
    assert M_int == pytest.approx(M, rel=1e-6)
    assert concrete_stress(sec, st).max() <= 0.0


def test_staged_design_is_the_smallest_strand_count_with_no_soffit_tension():
    d = design_prestress_staged(composite, IS_GIRDER, H_G, 98.7e-6, 1.0e9, M_A, M_A + DM)
    assert d.soffit_stress_service <= 0.0
    fewer = composite((d.n_strands - 1) * 98.7e-6)
    from crackfem.prestress import eps_p0_for_stress
    e0 = eps_p0_for_stress(girder_only(fewer, IS_GIRDER, H_G), 1.0e9)
    sec = staged_section(fewer, IS_GIRDER, H_G, e0, M_A)
    assert soffit_stress(sec, solve_state(sec, e0, M_A + DM, cracked=False)) > 0.0
