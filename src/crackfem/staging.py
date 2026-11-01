"""Staged construction of a precast prestressed girder with a cast-in-place slab (unshored).

Stage A: the bare girder carries the prestress and the weight of the girder and the wet slab.
Stage B: the hardened slab acts compositely; later loads act on the composite section.
The slab is cast stress-free on the deformed girder, so its stress-free strain is the
stage-A strain plane at its level. A composite solve at M = M_A therefore reproduces
stage A with an unstressed slab, and moments above M_A are stage-B increments.
"""
from dataclasses import replace
from typing import Callable

import numpy as np

from .geometry import Fibers
from .prestress import PrestressDesign, eps_p0_for_stress
from .section import PCSection, soffit_stress, solve_state


def girder_only(sec: PCSection, is_girder: np.ndarray, h_girder: float) -> PCSection:
    """The bare girder: the girder fibres of a composite section, same strand."""
    f = sec.fib
    fib = Fibers(f.y[is_girder], f.A[is_girder], f.E[is_girder], f.fct[is_girder])
    return PCSection(fib, sec.y_p, sec.A_p, sec.E_p, h_girder)


def staged_section(sec: PCSection, is_girder: np.ndarray, h_girder: float,
                   eps_p0_construction: float, M_A: float) -> PCSection:
    """Composite section whose slab fibres carry the stage-A strain plane as stress-free strain."""
    plane = solve_state(girder_only(sec, is_girder, h_girder), eps_p0_construction, M_A, cracked=False)
    return replace(sec, eps_free=np.where(is_girder, 0.0, plane.strain(sec.fib.y)))


def prestress_staged(make_section: Callable[[float], PCSection], is_girder: np.ndarray,
                     h_girder: float, n_strands: int, strand_area: float, fpe: float, M_A: float,
                     M_total: float) -> PrestressDesign:
    """A girder with a given strand count (as built): locked-in strain for fpe on the bare girder
    and, for the record, the staged soffit stress under the total service moment (> 0 = tension)."""
    if n_strands < 1:
        raise ValueError("n_strands must be >= 1")
    comp = make_section(n_strands * strand_area)
    e0 = eps_p0_for_stress(girder_only(comp, is_girder, h_girder), fpe)
    sec = staged_section(comp, is_girder, h_girder, e0, M_A)
    s = soffit_stress(sec, solve_state(sec, e0, M_total, cracked=False))
    return PrestressDesign(n_strands, n_strands * strand_area, fpe, e0, fpe * n_strands * strand_area,
                           M_total, s)


def design_prestress_staged(make_section: Callable[[float], PCSection], is_girder: np.ndarray,
                            h_girder: float, strand_area: float, fpe: float, M_A: float,
                            M_total: float, n_max: int = 150) -> PrestressDesign:
    """Smallest strand count with no soffit tension under the total service moment, with stage A
    on the bare girder. fpe = effective strand stress on the bare girder at zero external moment."""
    for n in range(1, n_max + 1):
        d = prestress_staged(make_section, is_girder, h_girder, n, strand_area, fpe, M_A, M_total)
        if d.soffit_stress_service <= 0.0:
            return d
    raise ValueError(f"no strand count up to {n_max} removes soffit tension")
