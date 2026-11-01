"""Direct stiffness assembly and partition solve of a grid (course method, Logan §5.4).

Member end forces are local: f^ = k^ T d - f^_eq = [V1, T1, M1, V2, T2, M2] acting on the member
(Logan's f^1y, m^1x, m^1z, f^2y, m^2x, m^2z). The internal sagging moment follows the beam rule
of solve.moment_at, M(x) = -M1 + V1 x + q x^2/2 + sum P (x - a), and the internal torque is -T1
(Logan Eq. 5.4.9)."""
from dataclasses import dataclass

import numpy as np
from numpy.linalg import LinAlgError
from scipy.linalg import cho_factor, cho_solve

from .elements import grid_k_local, grid_point_local, grid_T, grid_udl_local
from .gridmodel import GridModel


@dataclass(frozen=True)
class GridStiffness:
    K: np.ndarray            # assembled global stiffness
    kT: np.ndarray           # (n_members, 6, 6): k^ T of each member (local end forces = kT d - f^_eq)
    free: np.ndarray
    fixed: np.ndarray
    cho: tuple               # Cholesky factor of K_ff (K_ff is symmetric positive definite)


@dataclass(frozen=True)
class GridSolution:
    D: np.ndarray            # nodal displacements (w, theta_x, theta_z per node)
    R: np.ndarray            # reactions at restrained DOFs, zero elsewhere
    end_forces: np.ndarray   # (n_members, 6): [V1, T1, M1, V2, T2, M2], local


def member_dofs(model: GridModel) -> np.ndarray:
    i, j = model.members[:, 0], model.members[:, 1]
    return np.stack([3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2], axis=1)


def _per_member(values, model: GridModel, name: str) -> np.ndarray:
    v = np.broadcast_to(np.asarray(values, dtype=float), (model.n_members,))
    if not np.all(np.isfinite(v)) or np.any(v <= 0.0):
        raise ValueError(f"{name} must be finite and positive")
    return v


def stiffness(model: GridModel, EI, GJ) -> GridStiffness:
    EI, GJ = _per_member(EI, model, "EI"), _per_member(GJ, model, "GJ")
    L = model.lengths
    C, S = model.cosines
    dofs = member_dofs(model)
    K = np.zeros((model.n_dof, model.n_dof))
    kT = np.empty((model.n_members, 6, 6))
    for m in range(model.n_members):
        T = grid_T(C[m], S[m])
        kT[m] = grid_k_local(EI[m], GJ[m], L[m]) @ T
        K[np.ix_(dofs[m], dofs[m])] += T.T @ kT[m]            # Logan Eq. 5.4.19: k = T^T k^ T
    fixed = np.array(model.fixed_dofs, dtype=int)
    free = np.setdiff1d(np.arange(model.n_dof), fixed)
    try:
        cho = cho_factor(K[np.ix_(free, free)])
    except LinAlgError as exc:
        raise ValueError("K_ff not positive definite — check supports") from exc
    return GridStiffness(K, kT, free, fixed, cho)


def local_loads(model: GridModel) -> np.ndarray:
    """(n_members, 6) consistent local loads of each member's UDL and point loads."""
    L = model.lengths
    q = model.udl()
    f = np.array([grid_udl_local(q[m], L[m]) for m in range(model.n_members)])
    for m, a, P in model.point_loads:
        f[m] += grid_point_local(P, a, L[m])
    return f


def load_vector(model: GridModel) -> tuple[np.ndarray, np.ndarray]:
    """Global load vector F (consistent member loads rotated to global, plus nodal loads) and
    the local member loads."""
    C, S = model.cosines
    fl = local_loads(model)
    dofs = member_dofs(model)
    F = np.zeros(model.n_dof)
    for m in range(model.n_members):
        F[dofs[m]] += grid_T(C[m], S[m]).T @ fl[m]
    for dof, v in model.nodal_loads:
        F[dof] += v
    return F, fl


def solve(model: GridModel, EI=None, GJ=None, stiff: GridStiffness | None = None) -> GridSolution:
    """Partition solve K_ff D_f = F_f. Give either EI and GJ, or a GridStiffness of this model
    (to reuse its factorisation), not both."""
    if stiff is None:
        if EI is None or GJ is None:
            raise ValueError("give EI and GJ, or stiff")
        st = stiffness(model, EI, GJ)
    else:
        if EI is not None or GJ is not None:
            raise ValueError("give EI and GJ, or stiff, not both")
        if stiff.K.shape[0] != model.n_dof or stiff.kT.shape[0] != model.n_members:
            raise ValueError("stiff does not belong to this model (DOF or member count differs)")
        st = stiff
    F, fl = load_vector(model)
    D = np.zeros(model.n_dof)
    D[st.free] = cho_solve(st.cho, F[st.free])
    R = np.zeros(model.n_dof)
    R[st.fixed] = st.K[st.fixed] @ D - F[st.fixed]
    ef = np.einsum("mij,mj->mi", st.kT, D[member_dofs(model)]) - fl
    return GridSolution(D, R, ef)


def moment_at(model: GridModel, sol: GridSolution, m: int, xl: float) -> float:
    """Internal sagging moment at distance xl from node i of member m."""
    V1, M1 = sol.end_forces[m, 0], sol.end_forces[m, 2]
    out = -M1 + V1 * xl + model.udl()[m] * xl * xl / 2.0
    for mm, a, P in model.point_loads:
        if mm == m and a < xl:
            out += P * (xl - a)
    return float(out)


def torque(sol: GridSolution, m: int) -> float:
    """Internal torque of member m (constant along it), Logan Eq. 5.4.9."""
    return float(-sol.end_forces[m, 1])


def midpoint_moments(model: GridModel, sol: GridSolution, members) -> np.ndarray:
    L = model.lengths
    return np.array([moment_at(model, sol, int(m), L[int(m)] / 2.0) for m in members])


def line_node_moments(model: GridModel, sol: GridSolution, members) -> np.ndarray:
    """Sagging moment at every node of a chain of members (member k runs from node k to k + 1).
    The chain must be ordered and run in one direction: the end node of each member is the
    start node of the next (ValueError otherwise).
    At an interior node the moment jumps by the moment the cross members deliver there, so the
    node value is the mean of the two sides; the end nodes take their one side."""
    members = [int(m) for m in members]
    for a, b in zip(members[:-1], members[1:]):
        if model.members[a, 1] != model.members[b, 0]:
            raise ValueError(f"members {a} and {b} are not a same-direction chain "
                             f"(member {a} ends at node {model.members[a, 1]}, member {b} starts at node {model.members[b, 0]})")
    L = model.lengths
    left = [moment_at(model, sol, m, 0.0) for m in members]
    right = [moment_at(model, sol, m, L[m]) for m in members]
    inner = [(r + l) / 2.0 for r, l in zip(right[:-1], left[1:])]
    return np.array([left[0], *inner, right[-1]])
