import pytest

from crackfem.geometry import discretize, gross_properties, slab
from crackfem.section import PCSection, cracking_moment, resultants, solve_state, strand_stress, uncracked_EI

E_C, E_P = 30e9, 195e9
B, H, Y_P, A_P = 0.4, 0.8, 0.1, 1.5e-3


def rect_section(dy=2e-4):
    part = slab(0.0, H, B, E=E_C, fct=3e6)  # a rectangle is a trapezoid with b0 == b1
    return PCSection(discretize([part], dy), Y_P, A_P, E_P, H), part


def transformed():
    """Concrete rectangle + strand transformed by n = Ep/Ec (classic hand method)."""
    _, part = rect_section()
    Ac, yc, Ic = gross_properties([part])
    n = E_P / E_C
    At = Ac + n * A_P
    yt = (Ac * yc + n * A_P * Y_P) / At
    It = Ic + Ac * (yc - yt) ** 2 + n * A_P * (Y_P - yt) ** 2
    return At, yt, It


@pytest.mark.parametrize("M", [0.0, 150e3, 600e3])
def test_uncracked_state_is_in_equilibrium(M):
    sec, _ = rect_section()
    st = solve_state(sec, 4e-3, M, cracked=False)
    N, M_int = resultants(sec, 4e-3, st)
    assert abs(N) < 1e-6 * E_P * A_P * 4e-3
    assert M_int == pytest.approx(M, abs=1e-6 * 600e3)


def test_strand_stress_at_zero_moment_matches_transformed_section():
    sec, _ = rect_section()
    eps_p0 = 4e-3
    At, yt, It = transformed()
    P_lock = E_P * A_P * eps_p0                 # force the locked-in strand strain applies (not Ch3's {f0})
    e = yt - Y_P
    eps_c_at_strand = -P_lock / (E_C * At) - P_lock * e * e / (E_C * It)
    expected = E_P * (eps_p0 + eps_c_at_strand)
    got = strand_stress(sec, eps_p0, solve_state(sec, eps_p0, 0.0, cracked=False))
    assert got == pytest.approx(expected, rel=1e-6)


def test_uncracked_EI_matches_transformed_section():
    sec, _ = rect_section()
    _, _, It = transformed()
    assert uncracked_EI(sec) == pytest.approx(E_C * It, rel=1e-6)


def test_curvature_is_linear_in_moment_with_slope_one_over_EI():
    sec, _ = rect_section()
    k0 = solve_state(sec, 4e-3, 0.0, cracked=False).kappa
    k1 = solve_state(sec, 4e-3, 300e3, cracked=False).kappa
    assert (k1 - k0) == pytest.approx(300e3 / uncracked_EI(sec), rel=1e-9)


def test_cracking_moment_matches_hand_formula():
    sec, _ = rect_section()
    eps_p0, fct = 4e-3, 3e6
    At, yt, It = transformed()
    P_lock = E_P * A_P * eps_p0                 # force the locked-in strand strain applies (not Ch3's {f0})
    e = yt - Y_P
    # soffit stress: -P_lock/At - P_lock e yt / It + M yt / It = fct
    Mcr = (fct + P_lock / At + P_lock * e * yt / It) * It / yt
    assert cracking_moment(sec, eps_p0, fct) == pytest.approx(Mcr, rel=1e-5)
