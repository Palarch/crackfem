"""Expected Outcome 3 of the proposal: damage confined to a zone around midspan (the proposal's
one damaged element) against the distributed model, where every element takes its stiffness
from its own moment."""
import numpy as np
import pytest

from crackfem.girder1d import Girder1D
from crackfem.solve import solve

ETA = 0.36


@pytest.fixture(scope="module")
def g():
    return Girder1D()


@pytest.fixture(scope="module")
def meshes():
    return {n: Girder1D(n_el=n) for n in (10, 20, 40)}


def test_damaged_zone_marks_the_elements_within_half_its_length_of_midspan(g):
    L, n = g.inp.span, g.n_el
    pair = g.damaged_zone(2 * L / n)
    assert pair.sum() == 2 and pair[n // 2 - 1] and pair[n // 2]
    fifth = g.damaged_zone(L / 5)
    assert fifth.sum() * L / n == pytest.approx(L / 5, rel=1e-12)


def test_damage_everywhere_reproduces_the_distributed_model(g):
    full = g.indicator(ETA, damaged=np.ones(g.n_el, dtype=bool))
    assert full["deflection_truck_mm"] == pytest.approx(g.indicator(ETA)["deflection_truck_mm"], rel=1e-12)


def test_no_damaged_element_gives_the_linear_uncracked_response(g):
    none = g.indicator(ETA, damaged=np.zeros(g.n_el, dtype=bool))
    lin = solve(g.model_service, g.EI_I)
    assert none["deflection_service_mm"] == pytest.approx(-lin.D[2 * g.mid_node] * 1e3, rel=1e-9)
    assert min(none["EI_ratio"]) == pytest.approx(1.0, rel=1e-9)


def test_mesh_tied_damage_vanishes_as_the_mesh_is_refined(meshes):
    d = [m.indicator(ETA, damaged=m.damaged_zone(2 * m.inp.span / n))["deflection_truck_mm"]
         for n, m in meshes.items()]
    assert d[0] > d[1] > d[2]


def test_a_fixed_physical_length_converges_with_the_mesh(meshes):
    tied = {n: m.indicator(ETA, damaged=m.damaged_zone(2 * m.inp.span / n))["deflection_truck_mm"]
            for n, m in meshes.items()}
    fixed = {n: m.indicator(ETA, damaged=m.damaged_zone(m.inp.span / 5))["deflection_truck_mm"]
             for n, m in meshes.items()}
    assert abs(fixed[20] - fixed[40]) < 0.1 * abs(tied[20] - tied[40])
