"""Grid (grillage) model data (Logan 5th ed. §5.4). Nodes lie in the x-z plane, y is up.
DOFs per node (w, theta_x, theta_z) -> global (3i, 3i+1, 3i+2): w = vertical deflection (up +),
theta_x / theta_z = rotations about the global x / z axes (right-hand rule). Units N, m.
A member runs from node i to node j; its local x^ axis points from i to j."""
from dataclasses import dataclass, field, replace

import numpy as np


@dataclass(frozen=True)
class GridModel:
    xz: np.ndarray                                  # (n_nodes, 2): node coordinates (x, z)
    members: np.ndarray                             # (n_members, 2): node indices (i, j)
    fixed_dofs: tuple                               # restrained global DOFs
    q: np.ndarray | None = field(default=None)      # per-member uniform vertical load, up +
    point_loads: tuple = ()                         # ((member, a, P), ...): P up +, a from node i
    nodal_loads: tuple = ()                         # ((dof, value), ...): force or moment at a DOF

    def __post_init__(self):
        xz = np.asarray(self.xz, dtype=float)
        mem = np.asarray(self.members, dtype=int)
        if xz.ndim != 2 or xz.shape[1] != 2 or xz.shape[0] < 2:
            raise ValueError("xz must be an (n_nodes, 2) array with at least 2 nodes")
        if mem.ndim != 2 or mem.shape[1] != 2 or mem.shape[0] < 1:
            raise ValueError("members must be an (n_members, 2) array of node indices")
        if mem.min() < 0 or mem.max() >= xz.shape[0] or np.any(mem[:, 0] == mem[:, 1]):
            raise ValueError("member node indices must be distinct and in range")
        d = xz[mem[:, 1]] - xz[mem[:, 0]]
        if np.any(np.hypot(d[:, 0], d[:, 1]) <= 0.0):
            raise ValueError("members must have a positive length")
        n_dof = 3 * xz.shape[0]
        fixed = tuple(int(f) for f in self.fixed_dofs)
        if len(set(fixed)) != len(fixed) or any(not 0 <= f < n_dof for f in fixed):
            raise ValueError(f"fixed_dofs must be distinct, each in 0..{n_dof - 1}")
        object.__setattr__(self, "xz", xz)
        object.__setattr__(self, "members", mem)
        object.__setattr__(self, "fixed_dofs", fixed)
        # loads: the same rules whether set here or through with_loads
        if self.q is not None:
            try:
                q = np.broadcast_to(np.asarray(self.q, dtype=float), (mem.shape[0],)).copy()
            except ValueError as exc:
                raise ValueError("q must be a scalar or one value per member") from exc
            object.__setattr__(self, "q", q)
        L = self.lengths
        pts = tuple((int(m), float(a), float(P)) for m, a, P in self.point_loads)
        for m, a, _ in pts:
            if not 0 <= m < mem.shape[0] or not 0.0 <= a <= L[m]:
                raise ValueError(f"point load at a = {a} is not on member {m}")
        nod = tuple((int(dof), float(v)) for dof, v in self.nodal_loads)
        for dof, _ in nod:
            if not 0 <= dof < n_dof:
                raise ValueError(f"nodal load DOF {dof} is out of range")
        object.__setattr__(self, "point_loads", pts)
        object.__setattr__(self, "nodal_loads", nod)

    @property
    def n_nodes(self) -> int:
        return self.xz.shape[0]

    @property
    def n_dof(self) -> int:
        return 3 * self.n_nodes

    @property
    def n_members(self) -> int:
        return self.members.shape[0]

    @property
    def lengths(self) -> np.ndarray:
        d = self.xz[self.members[:, 1]] - self.xz[self.members[:, 0]]
        return np.hypot(d[:, 0], d[:, 1])

    @property
    def cosines(self) -> tuple[np.ndarray, np.ndarray]:
        """(C, S) per member: C = (xj - xi) / L, S = (zj - zi) / L (Logan Eq. 5.4.18)."""
        d = self.xz[self.members[:, 1]] - self.xz[self.members[:, 0]]
        L = self.lengths
        return d[:, 0] / L, d[:, 1] / L

    def udl(self) -> np.ndarray:
        return np.zeros(self.n_members) if self.q is None else self.q


def with_loads(model: GridModel, q=0.0, point_loads=(), nodal_loads=()) -> GridModel:
    """Copy of the model with these loads; replaces ALL existing loads. q: scalar or per member.
    Normalised and checked in GridModel.__post_init__."""
    return replace(model, q=q, point_loads=tuple(point_loads), nodal_loads=tuple(nodal_loads))
