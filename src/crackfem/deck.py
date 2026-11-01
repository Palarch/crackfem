"""The S01 deck as a grillage (Hambly's rules, sources S22-S27): girder lines along x, slab strips
across at every girder node station, diaphragms at both ends and at midspan.

Plan view: x along the span (0 at the first support), z across, girder g on the line z = g S
(G1 is g = 0). Node (g, k) sits at (x_k, z_g) and has index g (n_el + 1) + k. Members: girder g,
element k runs from node (g, k) to (g, k + 1); the cross member at station k between girders
g and g + 1 runs from node (g, k) to node (g + 1, k). Every girder end node is supported
vertically; rotations are free.

Loads between girder lines are shared by the lever rule (statically exact for vertical force
and for moment about the x axis); a load outside an exterior girder (on the overhang) goes to
that girder as a vertical load plus a torque about x: m_x = (z_load - z_girder) * P for a
downward load of magnitude P (right-hand rule, y up)."""
from dataclasses import dataclass

import numpy as np

from .gridmodel import GridModel, with_loads
from .materials import ec2_ecm
from .params import FT, INCH, Inputs


@dataclass(frozen=True)
class Deck:
    inp: Inputs
    n_el: int
    x: np.ndarray                  # station coordinates along the span (n_el + 1)
    z: np.ndarray                  # girder line coordinates (n_girders)
    model: GridModel               # geometry and supports, no loads
    girder_members: np.ndarray     # (n_girders, n_el) member ids
    cross_members: np.ndarray      # (n_el + 1, n_girders - 1) member ids
    kind: tuple                    # per member: "girder", "slab" or "diaphragm"
    strip_width: np.ndarray        # per member: tributary width along x of a cross member (0 for girders)

    @property
    def n_girders(self) -> int:
        return self.z.size

    @property
    def mid(self) -> int:
        return self.n_el // 2

    def node(self, g: int, k: int) -> int:
        return g * (self.n_el + 1) + k

    @property
    def overhang(self) -> float:
        return self.inp.overhang_width_ft * FT

    @property
    def edges(self) -> tuple[float, float]:
        return float(self.z[0] - self.overhang), float(self.z[-1] + self.overhang)

    @property
    def curbs(self) -> tuple[float, float]:
        bw = self.inp.barrier_width_ft * FT
        lo, hi = self.edges
        return lo + bw, hi - bw

    def girder_widths(self) -> np.ndarray:
        """Slab width of each composite girder: S for interior girders, S/2 + overhang outside."""
        S = self.inp.girder_spacing
        w = np.full(self.n_girders, S)
        w[0] = w[-1] = S / 2.0 + self.overhang
        return w

    def locate(self, x: float) -> tuple[int, float]:
        """(element index along a girder, distance from its first node) of a point at x."""
        L = self.inp.span
        tol = 1e-9 * L
        if x < -tol or x > L + tol:
            raise ValueError(f"x = {x} is outside the span")
        x = min(max(x, 0.0), L)
        k = min(int(np.searchsorted(self.x, x, side="right")) - 1, self.n_el - 1)
        return k, float(x - self.x[k])


def layout(inp: Inputs, n_el: int) -> Deck:
    if n_el < 2 or n_el % 2:
        raise ValueError("n_el must be even and >= 2 so that midspan is a station")
    ng = inp.n_girders
    if ng < 2:
        raise ValueError("a grillage needs at least 2 girders")
    x = np.linspace(0.0, inp.span, n_el + 1)
    z = inp.girder_spacing * np.arange(ng)
    xz = np.array([(xk, zg) for zg in z for xk in x])
    node = lambda g, k: g * (n_el + 1) + k
    members, kind, width = [], [], []
    for g in range(ng):
        for k in range(n_el):
            members.append((node(g, k), node(g, k + 1)))
            kind.append("girder")
            width.append(0.0)
    diaphragm_stations = {0, n_el // 2, n_el}
    dx = inp.span / n_el
    for k in range(n_el + 1):
        b = dx / 2.0 if k in (0, n_el) else dx
        for g in range(ng - 1):
            members.append((node(g, k), node(g + 1, k)))
            kind.append("diaphragm" if k in diaphragm_stations else "slab")
            width.append(b)
    fixed = tuple(3 * node(g, k) for g in range(ng) for k in (0, n_el))
    model = GridModel(xz, np.array(members), fixed)
    gm = np.arange(ng * n_el).reshape(ng, n_el)
    cm = ng * n_el + np.arange((n_el + 1) * (ng - 1)).reshape(n_el + 1, ng - 1)
    return Deck(inp, n_el, x, z, model, gm, cm, tuple(kind), np.array(width))


def shear_modulus(E: float, nu: float) -> float:
    return E / (2.0 * (1.0 + nu))


def torsion_constant_open(rects) -> float:
    """Logan Eq. 5.4.15: J = sum b t^3 / 3 over the thin rectangles (b, t) of an open section;
    t is taken as the smaller side of each rectangle."""
    return sum(max(b, t) * min(b, t) ** 3 / 3.0 for b, t in rects)


def girder_torsion_constant(inp: Inputs) -> float:
    """Type IV as three rectangles (Logan Table 5-1, wide-flange beam with unequal flanges):
    bottom flange, top flange, and a web between them; the haunches are left out."""
    parts = inp.girder_outline_in
    h_bot, b_bot = parts[0][0], parts[0][1]
    h_top, b_top = parts[-1][0], parts[-1][2]
    t_web = min(min(b0, b1) for _, b0, b1 in parts)
    h_web = sum(p[0] for p in parts) - h_bot - h_top
    return torsion_constant_open([(b_bot * INCH, h_bot * INCH), (b_top * INCH, h_top * INCH),
                                  (h_web * INCH, t_web * INCH)])


DIAPHRAGM_MODELS = ("web", "none", "composite")


def fixed_properties(deck: Deck, diaphragm: str = "web") -> tuple[np.ndarray, np.ndarray]:
    """(EI, GJ) per member. Cross members are linear: slab strip I = b t^3/12, J = b t^3/6
    (Hambly). Girder members get GJ = G_girder J_TypeIV + G_slab b t^3/6 with b = their slab
    width; their EI is NaN here (it comes from the section models).

    Diaphragms (S27 stand-in), three models for sensitivity:
    "web" (default): the web (depth x thickness, slab concrete; I = b h^3/12, J adds h b^3/3)
    added to the slab strip of its station, not composite with the slab, so its properties do not
    depend on the mesh. As a beam this is a lower bound on its stiffness; but the diaphragm is
    rigidly connected to the girder nodes (full moment transfer at the girder lines), which is an
    upper-bound assumption, so "web" is not a bound on the load sharing as a whole.
    "none": no intermediate (midspan) diaphragm - that station is a plain slab strip; the end
    diaphragms over the bearings stay as "web".
    "composite": every diaphragm is a T section, flange = the slab strip of its station, web below
    it (I about the T's own centroid; J as "web"). Its stiffness grows with the strip width L/n_el,
    so it depends on the mesh (the earlier prototype's model)."""
    if diaphragm not in DIAPHRAGM_MODELS:
        raise ValueError(f"diaphragm must be one of {DIAPHRAGM_MODELS}")
    inp = deck.inp
    t = inp.slab_thickness
    Es, Eg = ec2_ecm(inp.fck_slab), ec2_ecm(inp.fck_girder)
    Gs, Gg = shear_modulus(Es, inp.poisson_concrete), shear_modulus(Eg, inp.poisson_concrete)
    EI = np.full(deck.model.n_members, np.nan)
    GJ = np.empty(deck.model.n_members)
    Jg = girder_torsion_constant(inp)
    for g, w in enumerate(deck.girder_widths()):
        GJ[deck.girder_members[g]] = Gg * Jg + Gs * w * t**3 / 6.0
    hd, bd = inp.diaphragm_depth, inp.diaphragm_thickness
    mid_members = set(int(m) for m in deck.cross_members[deck.mid])
    for m in deck.cross_members.ravel():
        b = deck.strip_width[m]
        EI[m] = Es * b * t**3 / 12.0
        GJ[m] = Gs * b * t**3 / 6.0
        if deck.kind[m] != "diaphragm" or (diaphragm == "none" and int(m) in mid_members):
            continue
        GJ[m] += Gs * torsion_constant_open([(hd, bd)])
        if diaphragm == "composite":            # T: the strip as flange on top of the web
            A1, y1, A2, y2 = b * t, hd + t / 2.0, bd * hd, hd / 2.0
            yb = (A1 * y1 + A2 * y2) / (A1 + A2)
            EI[m] = Es * (b * t**3 / 12.0 + A1 * (y1 - yb) ** 2 + bd * hd**3 / 12.0 + A2 * (y2 - yb) ** 2)
        else:                                   # the web alone, added to the strip of its station
            EI[m] += Es * bd * hd**3 / 12.0
    return EI, GJ


def lever_weights(z_girders: np.ndarray, z: float) -> list[tuple[int, float, float]]:
    """[(girder, share of the load, torque arm)]: lever rule between two girder lines; outside
    the exterior girders the whole load with its arm z - z_girder."""
    zg = np.asarray(z_girders, dtype=float)
    if z <= zg[0]:
        return [(0, 1.0, float(z - zg[0]))]
    if z >= zg[-1]:
        return [(zg.size - 1, 1.0, float(z - zg[-1]))]
    i = int(np.searchsorted(zg, z, side="right")) - 1
    t = (z - zg[i]) / (zg[i + 1] - zg[i])
    return [(i, 1.0 - t, 0.0), (i + 1, t, 0.0)]


def line_to_girders(z_girders, z: float, p: float) -> tuple[np.ndarray, np.ndarray]:
    """A line load p (N/m, downward) along x at z -> (vertical line load per girder, downward;
    torque line load per girder, N·m/m about x)."""
    zg = np.asarray(z_girders, dtype=float)
    q, t = np.zeros(zg.size), np.zeros(zg.size)
    for g, share, arm in lever_weights(zg, z):
        q[g] += share * p
        t[g] += arm * share * p
    return q, t


def area_to_girders(z_girders, za: float, zb: float, w: float) -> tuple[np.ndarray, np.ndarray]:
    """A uniform area load w (N/m², downward) over za <= z <= zb -> per-girder line loads and
    torques. Exact: the strip is cut at the girder lines and each piece's resultant is shared
    by the lever rule at its centroid (the rule is linear in z between girder lines)."""
    zg = np.asarray(z_girders, dtype=float)
    cuts = np.unique(np.concatenate([[za, zb], zg[(zg > za) & (zg < zb)]]))
    q, t = np.zeros(zg.size), np.zeros(zg.size)
    for a, b in zip(cuts[:-1], cuts[1:]):
        dq, dt = line_to_girders(zg, 0.5 * (a + b), w * (b - a))
        q += dq
        t += dt
    return q, t


def stage_b_dead_lines(deck: Deck) -> tuple[np.ndarray, np.ndarray]:
    """Wearing surface between the curbs and the two barriers -> per-girder (line load, torque)."""
    inp = deck.inp
    lo, hi = deck.curbs
    q, t = area_to_girders(deck.z, lo, hi, inp.unit_weight_wearing * inp.wearing_thickness)
    bw = inp.barrier_width_ft * FT
    e_lo, e_hi = deck.edges
    for zb in (e_lo + bw / 2.0, e_hi - bw / 2.0)[: inp.n_barriers]:
        dq, dt = line_to_girders(deck.z, zb, inp.barrier_weight)
        q, t = q + dq, t + dt
    return q, t


def line_loads_model(deck: Deck, q_down: np.ndarray, t_line: np.ndarray, base: GridModel | None = None) -> GridModel:
    """Model with girder line loads (downward magnitudes) and their torques as consistent nodal
    moments about x (half of each element's share to each of its nodes)."""
    model = deck.model if base is None else base
    q = np.zeros(model.n_members)
    nodal = {}
    dx = deck.inp.span / deck.n_el
    for g in range(deck.n_girders):
        q[deck.girder_members[g]] = -q_down[g]
        if t_line[g] != 0.0:
            for k in range(deck.n_el + 1):
                share = dx / 2.0 if k in (0, deck.n_el) else dx
                dof = 3 * deck.node(g, k) + 1
                nodal[dof] = nodal.get(dof, 0.0) + t_line[g] * share
    return with_loads(model, q=q, nodal_loads=sorted(nodal.items()))


def wheel_loads(deck: Deck, axles, wheel_lines, scale: float) -> tuple[list, list]:
    """Point and nodal loads of a truck set: each axle (x, P) puts P/2 x scale on every wheel line
    z in wheel_lines (two lines per truck). Returns (point_loads, nodal_loads) in the GridModel
    form, up + (so the wheel loads are negative)."""
    points, nodal = [], {}
    for x, P in axles:
        k, a = deck.locate(float(x))
        L = deck.x[k + 1] - deck.x[k]
        for zw in wheel_lines:
            for g, share, arm in lever_weights(deck.z, float(zw)):
                load = share * P / 2.0 * scale            # downward magnitude
                points.append((int(deck.girder_members[g, k]), a, -load))
                if arm != 0.0:
                    for node, f in ((deck.node(g, k), 1.0 - a / L), (deck.node(g, k + 1), a / L)):
                        dof = 3 * node + 1
                        nodal[dof] = nodal.get(dof, 0.0) + arm * load * f
    return points, sorted(nodal.items())


def n_lanes(deck: Deck) -> int:
    lo, hi = deck.curbs
    return int(np.floor((hi - lo) / (deck.inp.lane_width_ft * FT) + 1e-9))


LANE_RULES = ("lanes", "clearance")


def truck_placements(deck: Deck, step: float = 0.05, lane_rule: str = "lanes") -> list[tuple[float, ...]]:
    """Transverse truck placements as tuples of wheel-line z: one truck, and (if the roadway has
    two or more design lanes) two trucks side by side. A truck occupies its clearance width
    (10 ft) with its wheel lines 2 ft inside it, and the clearance width stays between the curbs.

    lane_rule="lanes" (default, AASHTO Standard Specifications 3.6, S26): each of the two trucks
    sits inside its own 12 ft design lane; the two lanes do not overlap and lie anywhere between
    the curbs. lane_rule="clearance": the two 10 ft clearance widths only may not overlap (they
    may touch) - a closer, conservative pairing kept as a sensitivity. Edge positions are always
    included (curb positions; under "lanes" also a truck at the near edge of the second lane)."""
    if lane_rule not in LANE_RULES:
        raise ValueError(f"lane_rule must be one of {LANE_RULES}")
    inp = deck.inp
    lo, hi = deck.curbs
    cw = inp.truck_clearance_ft * FT
    edge = inp.wheel_curb_clearance_ft * FT
    sw = inp.wheel_spacing_ft * FT
    lane = inp.lane_width_ft * FT
    if hi - lo < cw:
        raise ValueError("the roadway is narrower than one truck")
    starts = np.unique(np.concatenate([np.arange(lo, hi - cw, step), [hi - cw]]))
    lines = lambda s: (s + edge, s + edge + sw)
    out = [lines(s) for s in starts]

    def pair_ok(s1: float, s2: float) -> bool:
        if s2 < s1 + cw - 1e-12:
            return False
        if lane_rule == "clearance":
            return True
        l1 = max(lo, s1 + cw - lane)                 # leftmost lane 1 that holds truck 1
        l2 = max(l1 + lane, s2 + cw - lane)          # leftmost lane 2 right of lane 1 that holds truck 2
        return l1 <= s1 + 1e-12 and l2 <= min(s2, hi - lane) + 1e-12

    if n_lanes(deck) >= 2 and hi - lo >= 2 * cw:
        pair_starts = starts
        if lane_rule == "lanes":                     # a truck at the near edge of the second lane
            pair_starts = np.unique(np.concatenate([starts, [lo + lane, hi - lane - cw]]))
        for s1 in pair_starts:
            out += [lines(s1) + lines(s2) for s2 in pair_starts if pair_ok(s1, s2)]
    # add the mirror image of every placement, so that a symmetric deck gets symmetric answers
    mirror = [tuple(sorted(lo + hi - z for z in p)) for p in out]
    rounded = {tuple(round(v, 9) for v in p): p for p in out + mirror}
    return [tuple(float(v) for v in rounded[key]) for key in sorted(rounded)]
