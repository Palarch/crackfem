"""Direct iteration with secant stiffness (course Ch11): solve, update each element's EI
from its own moment, solve again, until the moments stop changing."""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .model import BeamModel
from .solve import Solution, midpoint_moments, solve


@dataclass(frozen=True)
class NonlinearResult:
    solution: Solution
    EI: np.ndarray          # stiffness used in the final solve
    M_mid: np.ndarray       # element midpoint moments of the final solve
    iterations: int         # number of stiffness updates
    history: tuple          # relative moment change after each update


def solve_secant(model: BeamModel, stiffness_fn: Callable[[np.ndarray], np.ndarray], EI_start,
                 tol: float = 1e-6, max_iter: int = 100, tol_EI: float = 1e-6) -> NonlinearResult:
    """Stop when the moments have stopped changing (relative change <= tol) AND the stiffness
    belonging to the final moments equals the stiffness used to get them
    (max|stiffness_fn(M) - EI| / max(EI) <= tol_EI).

    stiffness_fn may be stateful (e.g. a sticky cracked flag). The stiffness-residual check
    prevents stopping with an EI that does not belong to the final moments, which a moment-only
    criterion allows when the stiffness jumps at Mcr. Raises RuntimeError if not converged."""
    EI = np.broadcast_to(np.asarray(EI_start, dtype=float), (model.n_el,)).copy()
    sol = solve(model, EI)
    M = midpoint_moments(model, sol)
    history = []
    for it in range(1, max_iter + 1):
        EI = np.asarray(stiffness_fn(M), dtype=float)
        sol = solve(model, EI)
        M_new = midpoint_moments(model, sol)
        change = float(np.max(np.abs(M_new - M)) / max(float(np.max(np.abs(M_new))), 1e-12))
        history.append(change)
        M = M_new
        if change <= tol:
            EI_check = np.asarray(stiffness_fn(M_new), dtype=float)
            if float(np.max(np.abs(EI_check - EI)) / np.max(EI)) <= tol_EI:
                return NonlinearResult(sol, EI, M, it, tuple(history))
    raise RuntimeError(f"secant iteration did not converge in {max_iter} updates: {history[-3:]}")


@dataclass(frozen=True)
class SecantResult:
    converged: bool
    solution: object        # whatever analyse returned for the final stiffness
    EI: np.ndarray          # stiffness used in the final solve
    M: np.ndarray           # moments of the final solve (what stiffness_fn reads)
    iterations: int         # number of stiffness updates
    history: tuple          # relative moment change after each update
    ei_residual: tuple      # max|stiffness_fn(M) - EI| / max(EI) after each update
    reason: str             # "converged" | "max_iter"


def secant_iterate(analyse: Callable[[np.ndarray], tuple], stiffness_fn: Callable[[np.ndarray], np.ndarray],
                   EI_start, tol: float = 1e-6, tol_EI: float = 1e-6, max_iter: int = 100,
                   relax: float = 1.0) -> SecantResult:
    """Direct iteration with secant stiffness for any FE model (Ch11): analyse(EI) returns
    (solution, M) and stiffness_fn(M) the secant EI belonging to M. Each update uses
    EI <- relax stiffness_fn(M) + (1 - relax) EI (relax = 1: plain secant; < 1: under-relaxed).
    Stops when the moments stop changing (relative change <= tol) and the stiffness belonging to
    the final moments equals the stiffness used (residual <= tol_EI). Never raises for
    non-convergence: returns converged=False with the whole history, so the caller can report it.
    stiffness_fn may keep state (e.g. a sticky pre-cracked flag) but must be idempotent for the
    same M: it is called twice per M (once for the residual check, once for the next update)."""
    if not 0.0 < relax <= 1.0:
        raise ValueError("relax must be in (0, 1]")
    EI = np.asarray(EI_start, dtype=float).copy()
    sol, M = analyse(EI)
    history, resid = [], []
    for it in range(1, max_iter + 1):
        EI = relax * np.asarray(stiffness_fn(M), dtype=float) + (1.0 - relax) * EI
        sol, M_new = analyse(EI)
        change = float(np.max(np.abs(M_new - M)) / max(float(np.max(np.abs(M_new))), 1e-12))
        r = float(np.max(np.abs(np.asarray(stiffness_fn(M_new), dtype=float) - EI)) / np.max(EI))
        history.append(change)
        resid.append(r)
        M = M_new
        if change <= tol and r <= tol_EI:
            return SecantResult(True, sol, EI, M, it, tuple(history), tuple(resid), "converged")
    return SecantResult(False, sol, EI, M, max_iter, tuple(history), tuple(resid), "max_iter")
