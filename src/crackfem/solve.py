"""Direct stiffness assembly and partition solve (course method)."""
from dataclasses import dataclass

import numpy as np

from .elements import beam_k, point_equivalent, udl_equivalent
from .model import BeamModel


@dataclass(frozen=True)
class Solution:
    D: np.ndarray             # nodal displacements (v, theta per node)
    R: np.ndarray             # reactions at restrained DOFs, zero elsewhere
    end_forces: np.ndarray    # (n_el, 4): [V1, M1, V2, M2] acting on each element


def _dofs(e: int) -> np.ndarray:
    return np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3])


def element_loads(model: BeamModel) -> list[np.ndarray]:
    q, L = model.udl(), model.lengths
    f = [udl_equivalent(q[e], L[e]) for e in range(model.n_el)]
    for xp, P in model.point_loads:
        e, a = model.locate(xp)
        f[e] = f[e] + point_equivalent(P, a, L[e])
    return f


def assemble(model: BeamModel, EI: np.ndarray):
    K = np.zeros((model.n_dof, model.n_dof))
    F = np.zeros(model.n_dof)
    fl = element_loads(model)
    for e in range(model.n_el):
        idx = _dofs(e)
        K[np.ix_(idx, idx)] += beam_k(EI[e], model.lengths[e])
        F[idx] += fl[e]
    return K, F, fl


def solve(model: BeamModel, EI) -> Solution:
    EI = np.broadcast_to(np.asarray(EI, dtype=float), (model.n_el,))
    if not np.all(np.isfinite(EI)) or np.any(EI <= 0):
        raise ValueError("EI must be finite and positive")
    K, F, fl = assemble(model, EI)
    fixed = np.array(model.fixed_dofs)
    free = np.setdiff1d(np.arange(model.n_dof), fixed)
    D = np.zeros(model.n_dof)
    D[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])   # partition: K_ff D_f = F_f
    R = np.zeros(model.n_dof)
    R[fixed] = K[fixed] @ D - F[fixed]
    ef = np.array([beam_k(EI[e], model.lengths[e]) @ D[_dofs(e)] - fl[e] for e in range(model.n_el)])
    return Solution(D, R, ef)


# Left free body of the element: the sagging moment acts counter-clockwise on the cut face, so M(0) = -M1 and M(L) = M2; summing moments about the cut gives M(x) = -M1 + V1 x + q x^2/2 + sum P (x - a) for loads with a < x.
def moment_at(model: BeamModel, sol: Solution, e: int, xl: float) -> float:
    """Internal sagging moment at distance xl from the first node of element e."""
    V1, M1 = sol.end_forces[e, 0], sol.end_forces[e, 1]
    m = -M1 + V1 * xl + model.udl()[e] * xl * xl / 2.0
    for xp, P in model.point_loads:
        ee, a = model.locate(xp)
        if ee == e and a < xl:
            m += P * (xl - a)
    return float(m)


def midpoint_moments(model: BeamModel, sol: Solution) -> np.ndarray:
    return np.array([moment_at(model, sol, e, model.lengths[e] / 2.0) for e in range(model.n_el)])


def node_moments(model: BeamModel, sol: Solution) -> np.ndarray:
    """Internal sagging moment at every node (start of each element, plus the last node)."""
    m = [moment_at(model, sol, e, 0.0) for e in range(model.n_el)]
    m.append(moment_at(model, sol, model.n_el - 1, model.lengths[-1]))
    return np.array(m)
