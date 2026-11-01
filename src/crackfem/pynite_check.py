"""Independent cross-check: the same linear grid solved by PyNite (PyNiteFEA, pinned in
requirements.txt). PyNite is a 3-D frame program, so the grid's in-plane DOFs (DX, DZ, RY) are
restrained at every node and the material is a unit one: E = G = 1, Iz = EI, J = GJ. For a
horizontal member PyNite's local y axis is the global Y and its local z = x cross y, as in
Logan's grid element, so PyNite's 'Mz' is the grid's bending moment (its sign convention is
PyNite's own; compare magnitudes or fit one sign)."""
import numpy as np

from .gridmodel import GridModel

COMBO = "C"
_DIRECTION = ("FY", "MX", "MZ")          # global DOF order (w, theta_x, theta_z)


def pynite_grid(model: GridModel, EI, GJ):
    """Build and solve the model in PyNite; returns the solved FEModel3D."""
    from Pynite.FEModel3D import FEModel3D

    EI = np.broadcast_to(np.asarray(EI, dtype=float), (model.n_members,))
    GJ = np.broadcast_to(np.asarray(GJ, dtype=float), (model.n_members,))
    fm = FEModel3D()
    for n, (x, z) in enumerate(model.xz):
        fm.add_node(f"N{n}", float(x), 0.0, float(z))
    fm.add_material("unit", 1.0, 1.0, 0.2, 0.0)
    for m, (i, j) in enumerate(model.members):
        fm.add_section(f"S{m}", 1.0, 1.0, float(EI[m]), float(GJ[m]))
        fm.add_member(f"M{m}", f"N{i}", f"N{j}", "unit", f"S{m}")
    fixed = set(model.fixed_dofs)
    for n in range(model.n_nodes):
        fm.def_support(f"N{n}", support_DX=True, support_DY=3 * n in fixed, support_DZ=True,
                       support_RX=3 * n + 1 in fixed, support_RY=True, support_RZ=3 * n + 2 in fixed)
    q = model.udl()
    for m in range(model.n_members):
        if q[m] != 0.0:
            fm.add_member_dist_load(f"M{m}", "FY", float(q[m]), float(q[m]))
    for m, a, P in model.point_loads:
        fm.add_member_pt_load(f"M{m}", "FY", float(P), float(a))
    for dof, v in model.nodal_loads:
        fm.add_node_load(f"N{dof // 3}", _DIRECTION[dof % 3], float(v))
    fm.add_load_combo(COMBO, {"Case 1": 1.0})
    fm.analyze_linear(check_statics=False)
    return fm


def pynite_deflections(fm, model: GridModel) -> np.ndarray:
    """Vertical deflection (up +) of every node."""
    return np.array([fm.nodes[f"N{n}"].DY[COMBO] for n in range(model.n_nodes)])


def pynite_moment(fm, member: int, x: float) -> float:
    """PyNite's bending moment 'Mz' of a member at distance x from its i node."""
    return float(fm.members[f"M{member}"].moment("Mz", x, COMBO))
