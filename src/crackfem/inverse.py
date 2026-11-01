"""Inverse problem: remaining prestress fraction eta from a measured crack width.

w_pred(eta) is the forward model (metres); it must not increase with eta. Widths at or
below w_tiny (1 nm) count as "no crack".
"""
import math
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .rootfind import illinois


@dataclass(frozen=True)
class EtaEstimate:
    status: str                    # "ok" | "no_crack_bound" | "exceeds_model"
    eta: float | None
    eta_lower_bound: float | None  # set for "no_crack_bound": eta >= this value
    iterations: int
    evaluations: int = 0           # number of w_pred calls this estimate cost


class _Counted:
    """Wraps w_pred: counts the calls and refuses non-finite widths."""

    def __init__(self, w_pred: Callable[[float], float]):
        self._w_pred = w_pred
        self.n = 0

    def __call__(self, eta: float) -> float:
        self.n += 1
        w = self._w_pred(eta)
        if not math.isfinite(w):
            raise ValueError(f"w_pred({eta!r}) is not finite")
        return w


def eta_decompression(w_pred: Callable[[float], float], eta_lo: float = 0.02, eta_hi: float = 1.0,
                      w_tiny: float = 1e-9, xtol: float = 1e-7) -> float:
    """Edge of the no-crack plateau: bisection on the predicate w_pred(eta) > w_tiny.

    Returns the cracked-side end of the final bracket, so w_pred(result) > w_tiny and the
    true edge lies within xtol above it. The predicate (not a root of w - w_tiny) is used
    because w_pred is flat at zero above the edge and may jump there."""
    w = w_pred if isinstance(w_pred, _Counted) else _Counted(w_pred)
    if w(eta_hi) > w_tiny:
        raise ValueError("the intact girder (eta_hi) already predicts a crack")
    if not (w(eta_lo) > w_tiny):
        return eta_lo
    lo, hi = eta_lo, eta_hi          # w(lo) > w_tiny (cracked), w(hi) <= w_tiny
    while hi - lo > xtol:
        mid = 0.5 * (lo + hi)
        if w(mid) > w_tiny:
            lo = mid
        else:
            hi = mid
    return lo


def estimate_eta(w_pred: Callable[[float], float], w_meas: float, eta_lo: float = 0.02,
                 eta_hi: float = 1.0, xtol: float = 1e-7, w_tiny: float = 1e-9) -> EtaEstimate:
    if not math.isfinite(w_meas):
        raise ValueError("w_meas is not finite")
    w = _Counted(w_pred)
    if w_meas <= w_tiny:
        bound = eta_decompression(w, eta_lo, eta_hi, w_tiny, xtol)
        return EtaEstimate("no_crack_bound", None, bound, 0, w.n)
    f = lambda e: w(e) - w_meas
    f_hi = f(eta_hi)
    if f_hi > 0.0:
        raise ValueError("the intact girder predicts a wider crack than measured: model inconsistent")
    f_lo = f(eta_lo)
    if f_lo < 0.0:
        return EtaEstimate("exceeds_model", None, None, 0, w.n)
    r = illinois(f, eta_lo, eta_hi, xtol=xtol, fa=f_lo, fb=f_hi)
    # If w_pred jumps (e.g. from a positive width straight to zero at decompression), the
    # bracket closes on the jump and f(r.x) is far from zero: no point estimate is honest.
    if abs(f(r.x)) > max(1e-3 * w_meas, w_tiny):
        return EtaEstimate("no_crack_bound", None, r.x, r.iterations, w.n)
    return EtaEstimate("ok", r.x, None, r.iterations, w.n)


def is_monotone_nonincreasing(w_pred: Callable[[float], float], grid,
                              tol: float | None = None) -> bool:
    grid = np.asarray(grid, dtype=float)
    if not np.all(np.diff(grid) > 0.0):
        raise ValueError("grid must be strictly ascending")
    w = np.array([w_pred(float(e)) for e in grid])
    if tol is None:
        wmax = float(np.max(np.abs(w))) if w.size else 0.0
        tol = 1e-9 * wmax if wmax > 0.0 else 1e-15
    return bool(np.all(np.diff(w) <= tol))
