"""Prestress sizing for the stand-in girder. eta = 1 means 'as designed'.

Design rule (source S06 + spec 5.1): the smallest number of strands such that the
intact girder has no soffit tension under the full service moment.
"""
from dataclasses import dataclass
from typing import Callable

from .section import PCSection, soffit_stress, solve_state, strand_stress


@dataclass(frozen=True)
class PrestressDesign:
    n_strands: int
    A_p: float
    fpe: float                 # effective strand stress at zero external moment
    eps_p0: float              # locked-in strand strain for eta = 1
    P0: float                  # strand force at zero external moment
    M_service: float
    soffit_stress_service: float


def effective_prestress(fpj_ratio: float, fpu: float, lump_sum_loss: float) -> float:
    return fpj_ratio * fpu - lump_sum_loss


def eps_p0_for_stress(sec: PCSection, fpe: float) -> float:
    """Locked-in strain giving strand stress fpe at zero external moment (linear in eps_p0)."""
    unit = strand_stress(sec, 1.0, solve_state(sec, 1.0, 0.0, cracked=False))
    return fpe / unit


def design_prestress(make_section: Callable[[float], PCSection], strand_area: float,
                     fpe: float, M_service: float, n_max: int = 150) -> PrestressDesign:
    for n in range(1, n_max + 1):
        sec = make_section(n * strand_area)
        e0 = eps_p0_for_stress(sec, fpe)
        s = soffit_stress(sec, solve_state(sec, e0, M_service, cracked=False))
        if s <= 0.0:
            return PrestressDesign(n, n * strand_area, fpe, e0, fpe * n * strand_area, M_service, s)
    raise ValueError(f"no strand count up to {n_max} removes soffit tension")
