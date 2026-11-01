import numpy as np
import sympy as sp

from crackfem.elements import beam_k, hermite, point_equivalent, udl_equivalent

x, L, EI, q, P, a = sp.symbols("x L EI q P a", positive=True)
xi = x / L
N = sp.Matrix([1 - 3 * xi**2 + 2 * xi**3,
               L * (xi - 2 * xi**2 + xi**3),
               3 * xi**2 - 2 * xi**3,
               L * (-xi**2 + xi**3)])
CASES = [(2.0e9, 1.3), (5.5, 0.4)]


def test_shape_functions_interpolate_nodal_values_and_slopes():
    dN = N.diff(x)
    assert list(N.subs(x, 0)) == [1, 0, 0, 0]
    assert list(dN.subs(x, 0)) == [0, 1, 0, 0]
    assert list(N.subs(x, L)) == [0, 0, 1, 0]
    assert list(sp.simplify(dN.subs(x, L))) == [0, 0, 0, 1]


def test_stiffness_equals_bending_energy_integral():
    B = N.diff(x, 2)
    k_sym = sp.lambdify((EI, L), sp.integrate(EI * B * B.T, (x, 0, L)))
    for EIv, Lv in CASES:
        np.testing.assert_allclose(beam_k(EIv, Lv), np.array(k_sym(EIv, Lv), float), rtol=1e-12)


def test_udl_equivalent_equals_work_integral():
    f_sym = sp.lambdify((q, L), sp.integrate(q * N, (x, 0, L)))
    np.testing.assert_allclose(udl_equivalent(-12.5, 1.7), np.array(f_sym(-12.5, 1.7), float).ravel(), rtol=1e-12)


def test_point_equivalent_equals_shape_functions_at_the_load():
    f_sym = sp.lambdify((P, a, L), P * N.subs(x, a))
    np.testing.assert_allclose(point_equivalent(-3.0, 0.45, 1.2),
                               np.array(f_sym(-3.0, 0.45, 1.2), float).ravel(), rtol=1e-12)
    np.testing.assert_allclose(hermite(0.45 / 1.2, 1.2),
                               np.array(sp.lambdify((a, L), N.subs(x, a))(0.45, 1.2), float).ravel(), rtol=1e-12)
