"""Tabulated forward map w(eta) for the Monte Carlo on the grillage.

Each exact w_pred(eta) is a nonlinear grillage solve, so the Monte Carlo (thousands of inversions)
runs on a table instead: w on a dense eta grid plus the edge of the no-crack plateau (found by
bisection on the exact map), joined by straight lines. Linear interpolation of non-increasing
data is non-increasing, so the inverse stays well posed; above the plateau edge the table returns
exactly zero, as the exact map does. The table is checked against the exact map (max_rel_error)
and the Monte Carlo percentiles are confirmed with the exact inverse (run_stage2). Below the
first grid point the table is flat (no extrapolation); the inverse's own bracket starts there too."""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .inverse import eta_decompression


@dataclass(frozen=True)
class WTable:
    eta: np.ndarray        # strictly ascending knots
    w: np.ndarray          # exact w_pred at the knots (m)
    eta_dec: float         # cracked-side edge of the no-crack plateau

    def __call__(self, eta: float) -> float:
        """Interpolated width; exactly 0 above eta_dec. Below the first knot (grid[0]) the table is
        flat at w(grid[0]) (np.interp clamps): it does not extrapolate."""
        if eta > self.eta_dec:
            return 0.0
        return float(np.interp(eta, self.eta, self.w))

    def is_nonincreasing(self) -> bool:
        return bool(np.all(np.diff(self.w) <= 1e-9 * max(float(self.w.max()), 1e-15)))


def build_table(w_pred: Callable[[float], float], grid, xtol: float = 1e-7) -> WTable:
    """Evaluate w_pred on the grid (ascending, within [grid[0], grid[-1]]) and at the plateau
    edge. Raises ValueError if any width is not finite (e.g. a non-converged secant loop) or the
    widths increase with eta (then the inverse is ill posed and the report must say so)."""
    grid = np.asarray(grid, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.diff(grid) > 0.0):
        raise ValueError("grid must be strictly ascending with at least 2 points")
    edge = eta_decompression(w_pred, float(grid[0]), float(grid[-1]), xtol=xtol)
    knots = np.unique(np.concatenate([grid[grid < edge], [edge]]))
    w = np.array([w_pred(float(e)) for e in knots])
    if not np.all(np.isfinite(w)):
        raise ValueError(f"w_pred is not finite at eta = {knots[~np.isfinite(w)].tolist()}")
    if knots[-1] < grid[-1]:                      # plateau: exact zeros from just above the edge
        knots = np.concatenate([knots, [min(edge + xtol, grid[-1])]])
        w = np.concatenate([w, [0.0]])
    t = WTable(knots, w, float(edge))
    if not t.is_nonincreasing():
        raise ValueError("w_pred increases with eta on the grid: the inverse is ill posed")
    return t


def max_rel_error(table: WTable, w_pred: Callable[[float], float], etas) -> float:
    """Largest |table - exact| / max(exact widths) over test points (e.g. the grid midpoints)."""
    exact = np.array([w_pred(float(e)) for e in etas])
    approx = np.array([table(float(e)) for e in etas])
    return float(np.max(np.abs(approx - exact)) / max(float(np.max(np.abs(exact))), 1e-15))
