"""Grid element (Logan 5th ed. §5.4): symbolic derivation and the printed Example 5.5."""
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from crackfem.elements import (beam_k, grid_k_global, grid_k_local, grid_point_local, grid_T,
                               grid_udl_local)
from crackfem.params import FT, INCH, KIP, KSI

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))
x, L, EI, GJ, q, P, a = sp.symbols("x L EI GJ q P a", positive=True)
xi = x / L
NW = sp.Matrix([1 - 3 * xi**2 + 2 * xi**3, L * (xi - 2 * xi**2 + xi**3), 3 * xi**2 - 2 * xi**3, L * (-xi**2 + xi**3)])
NT = sp.Matrix([1 - xi, xi])                         # linear twist, Logan Eq. 5.4.4


def _symbolic_k():
    kb = sp.integrate(EI * NW.diff(x, 2) * NW.diff(x, 2).T, (x, 0, L))
    kt = sp.integrate(GJ * NT.diff(x) * NT.diff(x).T, (x, 0, L))
    k = sp.zeros(6, 6)
    for r, R in enumerate((0, 2, 3, 5)):
        for c, Cc in enumerate((0, 2, 3, 5)):
            k[R, Cc] = kb[r, c]
    for r, R in enumerate((1, 4)):
        for c, Cc in enumerate((1, 4)):
            k[R, Cc] = kt[r, c]
    return sp.lambdify((EI, GJ, L), k)


def test_local_stiffness_equals_bending_plus_torsion_energy():
    k_sym = _symbolic_k()
    for EIv, GJv, Lv in ((2.0e9, 3.0e8, 1.3), (5.5, 0.7, 0.4)):
        np.testing.assert_allclose(grid_k_local(EIv, GJv, Lv), np.array(k_sym(EIv, GJv, Lv), float), rtol=1e-12)


def test_bending_block_is_the_beam_element():
    k = grid_k_local(3.0e9, 1.0e8, 2.0)
    idx = (0, 2, 3, 5)
    np.testing.assert_allclose(k[np.ix_(idx, idx)], beam_k(3.0e9, 2.0), rtol=1e-15)
    assert k[1, 1] == pytest.approx(1.0e8 / 2.0) and k[1, 4] == pytest.approx(-1.0e8 / 2.0)


def test_consistent_loads_are_the_work_integrals():
    fq = sp.lambdify((q, L), sp.integrate(q * NW, (x, 0, L)))
    fP = sp.lambdify((P, a, L), P * NW.subs(x, a))
    np.testing.assert_allclose(grid_udl_local(-12.5, 1.7)[[0, 2, 3, 5]], np.array(fq(-12.5, 1.7), float).ravel(), rtol=1e-12)
    np.testing.assert_allclose(grid_point_local(-3.0, 0.45, 1.2)[[0, 2, 3, 5]],
                               np.array(fP(-3.0, 0.45, 1.2), float).ravel(), rtol=1e-12)
    assert grid_udl_local(-12.5, 1.7)[1] == 0.0 and grid_point_local(-3.0, 0.45, 1.2)[4] == 0.0


@pytest.mark.parametrize("theta", [0.0, 0.3, 1.2, np.pi / 2, 2.5])
def test_transformation_is_orthogonal_and_keeps_the_deflection(theta):
    T = grid_T(np.cos(theta), np.sin(theta))
    np.testing.assert_allclose(T @ T.T, np.eye(6), atol=1e-15)
    assert T[0, 0] == 1.0 and T[3, 3] == 1.0


def test_global_stiffness_is_symmetric_with_rigid_body_modes():
    k = grid_k_global(2.0e9, 3.0e8, 2.5, 0.6, 0.8)
    np.testing.assert_allclose(k, k.T, rtol=1e-12, atol=1e-3)
    rigid = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])       # vertical translation
    np.testing.assert_allclose(k @ rigid, 0.0, atol=1e-3)


def _logan_example():
    ex = FIX["logan_grid_example"]
    E, G = ex["E_ksi"] * KSI, ex["G_ksi"] * KSI
    I, J = ex["I_in4"] * INCH**4, ex["J_in4"] * INCH**4
    xz = {int(n): np.array(c) * FT for n, c in ex["nodes_ft"].items()}
    return ex, E * I, G * J, xz


def test_logan_example_5_5_stiffness_and_displacements():
    ex, EIv, GJv, xz = _logan_example()
    K = np.zeros((3, 3))
    for i, j in ex["elements"]:
        d = xz[j] - xz[i]
        Lm = float(np.hypot(*d))
        K += grid_k_global(EIv, GJv, Lm, d[0] / Lm, d[1] / Lm)[:3, :3]
    # K in kip/in, kip, kip-in: rows/cols (w, phi_x, phi_z); w in in, phi in rad
    to_us = np.array([[KIP / INCH, KIP, KIP], [KIP, KIP * INCH, KIP * INCH], [KIP, KIP * INCH, KIP * INCH]])
    np.testing.assert_allclose(K / to_us, ex["K_free_kip_in"], rtol=5e-3, atol=1e-6)
    d = np.linalg.solve(K, [ex["load"]["F_y_kip"] * KIP, 0.0, 0.0])
    p = ex["displacements"]
    np.testing.assert_allclose([d[0] / INCH, d[1], d[2]], [p["d1y_in"], p["phi1x_rad"], p["phi1z_rad"]], rtol=5e-3)


def test_logan_example_5_5_local_element_forces():
    ex, EIv, GJv, xz = _logan_example()
    K = np.zeros((3, 3))
    geo = []
    for i, j in ex["elements"]:
        dd = xz[j] - xz[i]
        Lm = float(np.hypot(*dd))
        geo.append((Lm, dd[0] / Lm, dd[1] / Lm))
        K += grid_k_global(EIv, GJv, *geo[-1])[:3, :3]
    d1 = np.linalg.solve(K, [ex["load"]["F_y_kip"] * KIP, 0.0, 0.0])
    to_us = np.array([KIP, KIP * INCH, KIP * INCH, KIP, KIP * INCH, KIP * INCH])
    for e, (Lm, C, S) in enumerate(geo, start=1):
        f = grid_k_local(EIv, GJv, Lm) @ grid_T(C, S) @ np.concatenate([d1, np.zeros(3)])
        np.testing.assert_allclose(f / to_us, ex["local_forces_kip_kin"][str(e)], rtol=5e-3)
