import numpy as np
import pytest

from crackfem.geometry import discretize, slab
from crackfem.section import PCSection, cracking_moment, solve_state, uncracked_EI
from crackfem.stiffness import branson_secant_EI, ec2_secant_EI

SEC = PCSection(discretize([slab(0.0, 0.8, 0.4, E=30e9, fct=3e6)], 5e-4), 0.1, 1.5e-3, 195e9, 0.8)
E0 = 4e-3
MCR = cracking_moment(SEC, E0, 3e6)


def kI(M):
    return solve_state(SEC, E0, M, False).kappa


def kII(M):
    return solve_state(SEC, E0, M, True).kappa


def test_below_cracking_the_secant_stiffness_is_the_uncracked_stiffness():
    for M_ref in (0.0, 0.5 * MCR):
        assert ec2_secant_EI(SEC, E0, 0.9 * MCR, MCR, 0.5, M_ref) == pytest.approx(uncracked_EI(SEC), rel=1e-9)
        assert branson_secant_EI(SEC, E0, 0.9 * MCR, MCR, M_ref) == pytest.approx(uncracked_EI(SEC), rel=1e-9)


def test_ec2_interpolation_with_beta_one_at_twice_mcr():
    M = 2 * MCR
    expected = M / (0.75 * kII(M) + 0.25 * kI(M) - kI(0.0))      # zeta = 1 - (1/2)^2
    assert ec2_secant_EI(SEC, E0, M, MCR, 1.0) == pytest.approx(expected, rel=1e-12)


def test_secant_from_a_reference_moment_is_increment_over_increment():
    M, M_ref = 2 * MCR, 0.5 * MCR
    expected = (M - M_ref) / (0.75 * kII(M) + 0.25 * kI(M) - kI(M_ref))
    assert ec2_secant_EI(SEC, E0, M, MCR, 1.0, M_ref) == pytest.approx(expected, rel=1e-12)


def test_no_load_increment_keeps_the_uncracked_stiffness():
    assert ec2_secant_EI(SEC, E0, 1.5 * MCR, MCR, 0.5, M_ref=1.5 * MCR) == pytest.approx(uncracked_EI(SEC))


def test_secant_stiffness_falls_as_moment_grows_above_cracking():
    Ms = np.linspace(1.05 * MCR, 2.5 * MCR, 8)
    for M_ref in (0.0, 0.5 * MCR):
        ei = [ec2_secant_EI(SEC, E0, m, MCR, 0.5, M_ref) for m in Ms]
        assert all(b < a for a, b in zip(ei, ei[1:]))
        assert ei[0] < uncracked_EI(SEC)


def test_branson_lies_between_cracked_and_uncracked():
    M = 2 * MCR
    EI_cr = M / (kII(M) - kI(0.0))
    assert EI_cr < branson_secant_EI(SEC, E0, M, MCR) < uncracked_EI(SEC)


def test_reference_between_decompression_and_cracking_uses_kappa_one():
    M, M_ref = 2 * MCR, 0.9 * MCR
    expected = (M - M_ref) / (0.75 * kII(M) + 0.25 * kI(M) - kI(M_ref))
    assert ec2_secant_EI(SEC, E0, M, MCR, 1.0, M_ref) == pytest.approx(expected, rel=1e-12)


def test_cracked_reference_takes_both_curvatures_from_the_ec2_curve():
    M, M_ref = 2 * MCR, 1.25 * MCR
    z = 1 - 0.5 * (1 / 2) ** 2
    z_ref = 1 - 0.5 * (1 / 1.25) ** 2
    k = z * kII(M) + (1 - z) * kI(M)
    k_ref = z_ref * kII(M_ref) + (1 - z_ref) * kI(M_ref)
    assert ec2_secant_EI(SEC, E0, M, MCR, 0.5, M_ref) == pytest.approx((M - M_ref) / (k - k_ref), rel=1e-12)


def test_branson_formula_with_a_reference_moment():
    def kb(Mx):
        if Mx <= MCR:
            return kI(Mx)
        k0 = kI(0.0)
        r = (MCR / Mx) ** 3
        return k0 + Mx / (r * uncracked_EI(SEC) + (1 - r) * Mx / (kII(Mx) - k0))

    M, M_ref = 2 * MCR, 0.5 * MCR
    assert branson_secant_EI(SEC, E0, M, MCR, M_ref) == pytest.approx((M - M_ref) / (kb(M) - kb(M_ref)), rel=1e-12)
