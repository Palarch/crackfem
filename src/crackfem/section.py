"""Elastic fibre section of a bonded prestressed member.

Conventions (shared by all of crackfem):
  * y = height above the soffit (m); tension strain/stress positive.
  * kappa > 0 is sagging; strain eps(y) = eps0 - kappa * y (eps0 = soffit strain).
  * strand stress sigma_p = Ep * (eps_p0 + eps(y_p)); eps_p0 = strand strain locked
    in when the surrounding concrete has zero strain.
  * M = external moment (sagging +). Prestress is internal: the axial resultant
    of concrete + strand is zero.
Modes: cracked=False -> concrete linear in tension and compression;
       cracked=True  -> concrete carries compression only.
Staged sections: each fibre may carry a stress-free strain eps_free (the slab's strain locked in at casting); concrete stress is E * (eps(y) - eps_free).
"""
from dataclasses import dataclass

import numpy as np

from .geometry import Fibers


@dataclass(frozen=True)
class PCSection:
    fib: Fibers
    y_p: float
    A_p: float
    E_p: float
    h: float
    eps_free: np.ndarray | None = None      # per-fibre stress-free strain; None = all zero

    @property
    def soffit_E(self) -> float:
        return float(self.fib.E[int(np.argmin(self.fib.y))])

    def free_strain(self) -> np.ndarray:
        return np.zeros(self.fib.y.size) if self.eps_free is None else self.eps_free


@dataclass(frozen=True)
class State:
    eps0: float
    kappa: float
    M: float
    cracked: bool

    def strain(self, y):
        return self.eps0 - self.kappa * y

    @property
    def neutral_axis(self) -> float:
        """Height above the soffit where the strain is zero (needs kappa != 0)."""
        return self.eps0 / self.kappa


def _solve_active(sec: PCSection, eps_p0: float, M: float, active: np.ndarray) -> tuple[float, float]:
    f = sec.fib
    ea = f.E[active] * f.A[active]
    y = f.y[active]
    ef = sec.free_strain()[active]
    epap = sec.E_p * sec.A_p
    EA = ea.sum() + epap
    ES = (ea * y).sum() + epap * sec.y_p
    EI = (ea * y * y).sum() + epap * sec.y_p**2
    J = np.array([[EA, -ES], [-ES, EI]])
    rhs = np.array([-epap * eps_p0 + (ea * ef).sum(),
                    M + epap * sec.y_p * eps_p0 - (ea * ef * y).sum()])
    eps0, kappa = np.linalg.solve(J, rhs)
    return float(eps0), float(kappa)


def solve_state(sec: PCSection, eps_p0: float, M: float, cracked: bool) -> State:
    """Strain plane (eps0, kappa) that carries external moment M with zero axial force.
    Cracked mode: active-set iteration (fibres in tension drop out until the set is stable)."""
    if eps_p0 == 0.0 and M == 0.0 and sec.eps_free is None:
        return State(0.0, 0.0, 0.0, cracked)
    active = np.ones(sec.fib.y.size, dtype=bool)
    eps0, kappa = _solve_active(sec, eps_p0, M, active)
    if not cracked:
        return State(eps0, kappa, M, False)
    seen = set()
    for _ in range(200):
        new_active = (eps0 - kappa * sec.fib.y - sec.free_strain()) < 0.0
        if np.array_equal(new_active, active):
            return State(eps0, kappa, M, True)
        key = new_active.tobytes()
        if key in seen:
            raise RuntimeError("cracked-section active set cycled")
        seen.add(key)
        if not new_active.any():
            raise ValueError("no concrete left in compression")
        active = new_active
        eps0, kappa = _solve_active(sec, eps_p0, M, active)
    raise RuntimeError("cracked-section active set did not converge")


def concrete_stress(sec: PCSection, st: State) -> np.ndarray:
    sig = sec.fib.E * (st.strain(sec.fib.y) - sec.free_strain())
    return np.minimum(sig, 0.0) if st.cracked else sig


def strand_stress(sec: PCSection, eps_p0: float, st: State) -> float:
    return float(sec.E_p * (eps_p0 + st.strain(sec.y_p)))


def delta_sigma_p(sec: PCSection, st: State) -> float:
    """Strand stress change from the state of zero concrete strain at the strand level."""
    return float(sec.E_p * st.strain(sec.y_p))


def resultants(sec: PCSection, eps_p0: float, st: State) -> tuple[float, float]:
    """(axial force, internal sagging moment) of the stresses in state st."""
    sig = concrete_stress(sec, st)
    sp = strand_stress(sec, eps_p0, st)
    N = float((sig * sec.fib.A).sum() + sp * sec.A_p)
    M = float(-((sig * sec.fib.A * sec.fib.y).sum() + sp * sec.A_p * sec.y_p))
    return N, M


def uncracked_EI(sec: PCSection) -> float:
    f = sec.fib
    ea = f.E * f.A
    epap = sec.E_p * sec.A_p
    EA = ea.sum() + epap
    ES = (ea * f.y).sum() + epap * sec.y_p
    EI = (ea * f.y * f.y).sum() + epap * sec.y_p**2
    return float(EI - ES * ES / EA)


def soffit_stress(sec: PCSection, st: State) -> float:
    """Concrete stress at y = 0 for an uncracked state."""
    if st.cracked:
        raise ValueError("soffit_stress is defined for uncracked states only")
    return sec.soffit_E * (st.eps0 - float(sec.free_strain()[int(np.argmin(sec.fib.y))]))


def cracking_moment(sec: PCSection, eps_p0: float, fct: float) -> float:
    """External moment at which the uncracked soffit stress reaches fct (fct = 0 gives the
    decompression moment). Soffit stress is linear in M, so two solves suffice."""
    s0 = soffit_stress(sec, solve_state(sec, eps_p0, 0.0, cracked=False))
    s1 = soffit_stress(sec, solve_state(sec, eps_p0, 1.0e6, cracked=False))
    return 1.0e6 * (fct - s0) / (s1 - s0)
