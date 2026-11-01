import pytest

from crackfem.geometry import discretize, slab
from crackfem.prestress import design_prestress, effective_prestress, eps_p0_for_stress
from crackfem.section import PCSection, solve_state, soffit_stress, strand_stress

FIB = discretize([slab(0.0, 0.8, 0.4, E=30e9, fct=3e6)], 5e-4)


def make_section(A_p):
    return PCSection(FIB, 0.1, A_p, 195e9, 0.8)


def test_effective_prestress_is_jacking_minus_lump_sum_loss():
    assert effective_prestress(0.75, 1860e6, 300e6) == pytest.approx(1095e6)


def test_locked_in_strain_reproduces_the_target_stress_at_zero_moment():
    sec = make_section(1.0e-3)
    e0 = eps_p0_for_stress(sec, 1.0e9)
    assert strand_stress(sec, e0, solve_state(sec, e0, 0.0, cracked=False)) == pytest.approx(1.0e9)


def test_design_is_the_smallest_strand_count_with_no_soffit_tension():
    d = design_prestress(make_section, 98.7e-6, 1.0e9, 400e3)
    assert d.soffit_stress_service <= 0.0
    fewer = make_section((d.n_strands - 1) * 98.7e-6)
    e0 = eps_p0_for_stress(fewer, 1.0e9)
    assert soffit_stress(fewer, solve_state(fewer, e0, 400e3, cracked=False)) > 0.0
    assert d.P0 == pytest.approx(d.n_strands * 98.7e-6 * 1.0e9)
