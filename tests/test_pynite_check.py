"""Independent solver: the same linear grid in PyNite (PyNiteFEA, pinned in requirements.txt)."""
from pathlib import Path

import numpy as np
import pytest

from crackfem.gridmodel import GridModel, with_loads
from crackfem.gridsolve import moment_at, solve
from crackfem.grillage import Grillage
from crackfem.pynite_check import pynite_deflections, pynite_grid, pynite_moment

REQ = Path(__file__).resolve().parents[1] / "requirements.txt"


def test_pynite_is_pinned():
    lines = [line.strip().lower() for line in REQ.read_text(encoding="utf-8").splitlines()]
    assert "pynitefea==2.0.2" in lines                         # pip may print the name in either case


def test_simply_supported_member_against_the_closed_form():
    L, EI, q = 8.0, 2.0e7, -10e3
    s = np.linspace(0.0, L, 3)
    model = with_loads(GridModel(np.stack([s, np.zeros(3)], axis=1), np.array([(0, 1), (1, 2)]), (0, 1, 6)), q=q)
    fm = pynite_grid(model, EI, 1.0e7)
    w = pynite_deflections(fm, model)
    assert w[1] == pytest.approx(5 * q * L**4 / (384 * EI), rel=1e-6)
    assert abs(pynite_moment(fm, 1, 0.0)) == pytest.approx(abs(q) * L**2 / 8, rel=1e-6)


@pytest.fixture(scope="module")
def G():
    return Grillage(n_el=10)


def test_intact_grillage_service_case_agrees_with_pynite(G):
    model = G.service_model(0)
    ours = solve(model, G.EI_I, G.GJ0)
    fm = pynite_grid(model, G.EI_I, G.GJ0)
    w = pynite_deflections(fm, model)
    np.testing.assert_allclose(w, ours.D[0::3], rtol=1e-6, atol=1e-6 * np.abs(ours.D[0::3]).max())
    d = G.deck
    right_of_mid = [int(d.girder_members[g, d.mid]) for g in range(d.n_girders)]   # starts at midspan
    ours_m = np.array([moment_at(model, ours, m, 0.0) for m in right_of_mid])
    pyn_m = np.array([pynite_moment(fm, m, 0.0) for m in right_of_mid])
    sign = np.sign(pyn_m[2] / ours_m[2])                         # PyNite's own Mz sign convention
    np.testing.assert_allclose(sign * pyn_m, ours_m, rtol=1e-6, atol=1e-6 * np.abs(ours_m).max())
