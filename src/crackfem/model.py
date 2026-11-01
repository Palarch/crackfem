"""1-D beam model data. DOFs per node (v, theta) -> global (2i, 2i+1). Units N, m."""
from dataclasses import dataclass, field, replace

import numpy as np


@dataclass(frozen=True)
class BeamModel:
    x: np.ndarray                                   # node coordinates
    fixed_dofs: tuple                               # restrained global DOFs
    q: np.ndarray | None = field(default=None)      # per-element uniform load, up +
    point_loads: tuple = ()                         # ((x, P), ...), up +

    def __post_init__(self):
        x = np.asarray(self.x, dtype=float)
        if x.ndim != 1 or x.size < 2 or not np.all(np.diff(x) > 0.0):
            raise ValueError("x must be a strictly increasing 1-D array of at least 2 nodes")
        object.__setattr__(self, "x", x)
        n_dof = 2 * x.size
        if len(self.fixed_dofs) < 2 or any(not 0 <= d < n_dof for d in self.fixed_dofs):
            raise ValueError(f"fixed_dofs needs at least 2 restraints, each in 0..{n_dof - 1}")

    @property
    def n_el(self) -> int:
        return self.x.size - 1

    @property
    def n_dof(self) -> int:
        return 2 * self.x.size

    @property
    def lengths(self) -> np.ndarray:
        return np.diff(self.x)

    @property
    def midpoints(self) -> np.ndarray:
        return 0.5 * (self.x[:-1] + self.x[1:])

    def udl(self) -> np.ndarray:
        return np.zeros(self.n_el) if self.q is None else self.q

    def locate(self, xp: float) -> tuple[int, float]:
        """(element index, distance from its first node) of a point on the beam."""
        tol = 1e-9 * (self.x[-1] - self.x[0])
        if xp < self.x[0] - tol or xp > self.x[-1] + tol:
            raise ValueError(f"x = {xp} is outside the beam")
        xp = min(max(xp, self.x[0]), self.x[-1])      # clip rounding overshoot to the end
        e = min(int(np.searchsorted(self.x, xp, side="right")) - 1, self.n_el - 1)
        return e, float(xp - self.x[e])


def simply_supported(L: float, n_el: int) -> BeamModel:
    return BeamModel(x=np.linspace(0.0, L, n_el + 1), fixed_dofs=(0, 2 * n_el))


def with_loads(model: BeamModel, q: float = 0.0, point_loads=()) -> BeamModel:
    """Return a copy of the model with these loads; replaces ALL existing loads (pass the UDL and the point loads in one call)."""
    return replace(model, q=np.full(model.n_el, float(q)), point_loads=tuple(point_loads))
