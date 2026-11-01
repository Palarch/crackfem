"""Euler-Bernoulli beam element (course Ch3). DOFs per element: (v1, theta1, v2, theta2);
v up positive, theta counter-clockwise positive, loads up positive."""
import numpy as np


def beam_k(EI: float, L: float) -> np.ndarray:
    c = EI / L**3
    return c * np.array([
        [12.0, 6.0 * L, -12.0, 6.0 * L],
        [6.0 * L, 4.0 * L * L, -6.0 * L, 2.0 * L * L],
        [-12.0, -6.0 * L, 12.0, -6.0 * L],
        [6.0 * L, 2.0 * L * L, -6.0 * L, 4.0 * L * L],
    ])


def hermite(xi: float, L: float) -> np.ndarray:
    """Cubic Hermite shape functions at xi = x / L."""
    return np.array([1 - 3 * xi**2 + 2 * xi**3,
                     L * (xi - 2 * xi**2 + xi**3),
                     3 * xi**2 - 2 * xi**3,
                     L * (-xi**2 + xi**3)])


def udl_equivalent(q: float, L: float) -> np.ndarray:
    """Consistent nodal loads of a uniform load q (up +) over the element."""
    return np.array([q * L / 2.0, q * L * L / 12.0, q * L / 2.0, -q * L * L / 12.0])


def point_equivalent(P: float, a: float, L: float) -> np.ndarray:
    """Consistent nodal loads of a point load P (up +) at distance a from node 1."""
    return P * hermite(a / L, L)


# --- Grid element (Logan 5th ed. §5.4, Eqs. 5.4.16-5.4.19) -------------------------------
# Local DOFs per element: (w1, phi1x, phi1z, w2, phi2x, phi2z): vertical deflection (y up +),
# twist about the member axis x^, bending rotation about z^ = x^ cross y^ (the beam's theta).
# The bending block is beam_k on (w, phi_z); the torsion block is GJ/L [[1, -1], [-1, 1]].
GRID_BENDING = (0, 2, 3, 5)
GRID_TORSION = (1, 4)


def grid_k_local(EI: float, GJ: float, L: float) -> np.ndarray:
    k = np.zeros((6, 6))
    k[np.ix_(GRID_BENDING, GRID_BENDING)] = beam_k(EI, L)
    t = GJ / L
    k[np.ix_(GRID_TORSION, GRID_TORSION)] = t * np.array([[1.0, -1.0], [-1.0, 1.0]])
    return k


def grid_T(C: float, S: float) -> np.ndarray:
    """Logan Eq. 5.4.18: local = T @ global per node, with C = (xj - xi) / L, S = (zj - zi) / L."""
    r = np.array([[1.0, 0.0, 0.0], [0.0, C, S], [0.0, -S, C]])
    T = np.zeros((6, 6))
    T[:3, :3] = r
    T[3:, 3:] = r
    return T


def grid_k_global(EI: float, GJ: float, L: float, C: float, S: float) -> np.ndarray:
    """Logan Eq. 5.4.19: k = T^T k^ T."""
    T = grid_T(C, S)
    return T.T @ grid_k_local(EI, GJ, L) @ T


def grid_udl_local(q: float, L: float) -> np.ndarray:
    """Consistent local loads of a uniform vertical load q (up +): the beam's, no torque."""
    f = np.zeros(6)
    f[list(GRID_BENDING)] = udl_equivalent(q, L)
    return f


def grid_point_local(P: float, a: float, L: float) -> np.ndarray:
    """Consistent local loads of a vertical point load P (up +) at distance a from node 1."""
    f = np.zeros(6)
    f[list(GRID_BENDING)] = point_equivalent(P, a, L)
    return f
