"""Section outline as stacked trapezoids, exact gross properties, fibre discretisation.

Elevation y is measured upward from the girder soffit (m).
"""
import math
from dataclasses import dataclass

import numpy as np

from .params import INCH


@dataclass(frozen=True)
class Trapezoid:
    y0: float          # bottom elevation
    h: float           # height
    b0: float          # width at y0
    b1: float          # width at y0 + h
    E: float           # modulus of this concrete
    fct: float         # tensile strength of this concrete
    is_girder: bool = True

    @property
    def area(self) -> float:
        return self.h * (self.b0 + self.b1) / 2.0

    @property
    def centroid(self) -> float:
        return self.y0 + self.h * (self.b0 + 2.0 * self.b1) / (3.0 * (self.b0 + self.b1))

    @property
    def inertia_own(self) -> float:
        b0, b1 = self.b0, self.b1
        return self.h**3 * (b0 * b0 + 4.0 * b0 * b1 + b1 * b1) / (36.0 * (b0 + b1))

    def width_at(self, y):
        return self.b0 + (self.b1 - self.b0) * (y - self.y0) / self.h

    def area_below(self, y_cut: float) -> float:
        if y_cut <= self.y0:
            return 0.0
        if y_cut >= self.y0 + self.h:
            return self.area
        t = y_cut - self.y0
        return t * (self.b0 + self.width_at(y_cut)) / 2.0


def girder_outline(dims_in, E: float, fct: float) -> list[Trapezoid]:
    """Stack (height_in, width_bottom_in, width_top_in) parts from the soffit upward."""
    parts, y = [], 0.0
    for h_in, b0_in, b1_in in dims_in:
        h = h_in * INCH
        parts.append(Trapezoid(y, h, b0_in * INCH, b1_in * INCH, E, fct, True))
        y += h
    return parts


def slab(y0: float, thickness: float, width: float, E: float, fct: float) -> Trapezoid:
    return Trapezoid(y0, thickness, width, width, E, fct, False)


def top_of(parts) -> float:
    return max(p.y0 + p.h for p in parts)


def gross_properties(parts, E_ref=None) -> tuple[float, float, float]:
    """(area, centroid height, second moment about the centroid); each part is
    weighted by E / E_ref when E_ref is given (transformed section)."""
    w = [1.0 if E_ref is None else p.E / E_ref for p in parts]
    A = sum(wi * p.area for wi, p in zip(w, parts))
    ybar = sum(wi * p.area * p.centroid for wi, p in zip(w, parts)) / A
    I = sum(wi * (p.inertia_own + p.area * (p.centroid - ybar) ** 2) for wi, p in zip(w, parts))
    return A, ybar, I


def area_below(parts, y_cut: float, girder_only: bool = True) -> float:
    """Plain (not E-weighted) area of the parts below y_cut; girder parts only by default.
    Used for EC2's effective tension area Ac,eff — always use this geometric value, not a fibre sum."""
    return sum(p.area_below(y_cut) for p in parts if p.is_girder or not girder_only)


@dataclass(frozen=True)
class Fibers:
    y: np.ndarray
    A: np.ndarray
    E: np.ndarray
    fct: np.ndarray


def discretize(parts, dy_max: float = 1e-3) -> Fibers:
    """Horizontal layers no thicker than dy_max; layer area = width at mid-height x thickness
    (exact for a linearly varying width)."""
    ys, As, Es, fs = [], [], [], []
    for p in parts:
        n = max(1, math.ceil(p.h / dy_max))
        dy = p.h / n
        y_mid = p.y0 + dy * (np.arange(n) + 0.5)
        ys.append(y_mid)
        As.append(p.width_at(y_mid) * dy)
        Es.append(np.full(n, p.E))
        fs.append(np.full(n, p.fct))
    return Fibers(np.concatenate(ys), np.concatenate(As), np.concatenate(Es), np.concatenate(fs))
