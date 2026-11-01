"""Stage 2: the chain on the grillage of the S01 deck (spec §3, Stage 2).

Each girder g is a Stage-1 Girder1D built with its own slab width (exterior: S/2 + overhang;
interior: S) and its own strand count, so its staged sections, stage-A moments (statics on the
bare girder), Mcr(eta) and crack width are exactly the Stage-1 physics (reused, not copied).
Stage B (wearing surface, barriers, truck) acts on the grillage. The secant loop updates the EI
of every girder element from that element's own total moment; cross members and torsion stay
linear. The cracked girder Gk has eta < 1, the others eta = 1 under the same rule.

Pre-cracked elements: an element whose total service moment under the moving-truck envelope
(linear intact grillage, the truck at girder g's critical transverse placement) exceeds Mcr(eta)
starts pre-cracked (fct = 0, as Stage 1).

Loading path: stage B is applied in increments - the stage-B dead load first, then the truck in
n_steps equal steps - and the secant loop converges at every step (warm-started from the
previous step). When a converged state puts unflagged elements above their sound Mcr, only the
one with the largest M/Mcr is flagged pre-cracked (decompression Mcr) for the rest of the path
(sticky flag) and the step is re-converged - one cracking event at a time, so the others are
re-checked after the shed load has been redistributed. Flags are set from converged equilibrium
states only, never from an intermediate iterate, so they do not depend on the iteration path
(relax), and one event at a time makes them independent of the load-step size (n_steps).

At-rest (w_pred) analysis: the stage-B dead load alone, starting from the envelope flags above.
Flags that a damaged service run would add to the neighbours are not carried into it: the truck
history is represented by the linear intact moving-truck envelope (the same rule for every
girder), and a single damaged service run is one truck position, not a history."""
from dataclasses import dataclass, replace

import numpy as np

from .deck import (fixed_properties, girder_torsion_constant, layout, lever_weights, line_loads_model,
                   stage_b_dead_lines, truck_placements, wheel_loads)
from .girder1d import Girder1D
from .gridmodel import with_loads
from .gridsolve import line_node_moments, midpoint_moments, solve, stiffness
from .loads import axle_offsets, impact_factor, truck_position_max_midspan
from .nonlinear import SecantResult, secant_iterate
from .params import Inputs
from .prestress import effective_prestress
from .staging import design_prestress_staged, prestress_staged
from .stiffness import ec2_secant_EI


@dataclass(frozen=True)
class GrillageAnalysis:
    result: SecantResult          # secant result of the last load step (full load), or of the step that failed
    etas: tuple
    precracked: np.ndarray        # (n_girders, n_el) after the analysis (sticky flags included)
    M_B_mid: np.ndarray           # stage-B sagging moment at each girder's midspan node
    deflection_mid: np.ndarray    # downward midspan deflection of each girder (m)
    step_iterations: tuple = ()   # secant updates per load step (re-solves after new flags included)
    # (n_girders, n_el) largest total moment along the path: over accepted converged states (no
    # unflagged element above its sound Mcr) and, for a flagged element, at its cracking event
    path_max: np.ndarray | None = None


def _same_loads(a, b) -> bool:
    return (np.array_equal(a.udl(), b.udl()) and tuple(a.point_loads) == tuple(b.point_loads)
            and tuple(a.nodal_loads) == tuple(b.nodal_loads))


def _blend(a, b, s: float):
    """Model a with the loads (1 - s) a + s b (loads are linear, so this is a + s (b - a))."""
    nodal = {}
    for w, m in ((1.0 - s, a), (s, b)):
        for dof, v in m.nodal_loads:
            nodal[dof] = nodal.get(dof, 0.0) + w * v
    pts = [(m, x, (1.0 - s) * P) for m, x, P in a.point_loads] + [(m, x, s * P) for m, x, P in b.point_loads]
    return with_loads(a, q=(1.0 - s) * a.udl() + s * b.udl(), point_loads=pts, nodal_loads=sorted(nodal.items()))


class Grillage:
    def __init__(self, inp: Inputs = Inputs(), n_el: int = 20, dy_max: float = 1e-3,
                 n_strands=None, placement_step: float = 0.05, max_design_rounds: int = 6,
                 n_steps: int = 10, lane_rule: str = "lanes", diaphragm: str = "web"):
        """n_strands=None designs every girder by the Stage-1 zero-tension rule under its own
        grillage service moment (eta = 1); a tuple analyses those as-built counts. n_steps: equal
        truck increments of the nonlinear load path (after the stage-B dead load). lane_rule:
        two-truck placement rule of `deck.truck_placements` ("lanes" = AASHTO 12 ft design lanes;
        "clearance" = touching 10 ft clearance widths, a conservative sensitivity). diaphragm:
        "web" | "none" | "composite", the diaphragm model of `deck.fixed_properties` (sensitivity)."""
        if n_steps < 1:
            raise ValueError("n_steps must be >= 1")
        self.inp, self.n_el, self.dy_max, self.n_steps = inp, n_el, dy_max, int(n_steps)
        self.lane_rule, self.diaphragm = lane_rule, diaphragm
        self.deck = layout(inp, n_el)
        self.widths = self.deck.girder_widths()
        self.EI_fixed, self.GJ0 = fixed_properties(self.deck, diaphragm)
        self.J_girder = girder_torsion_constant(inp)
        self.q_dead, self.t_dead = stage_b_dead_lines(self.deck)
        self.model_dead = line_loads_model(self.deck, self.q_dead, self.t_dead)
        self.axles, _ = truck_position_max_midspan(inp)         # Stage-1 longitudinal placement
        self.truck_scale = (1.0 + impact_factor(inp)) * inp.live_load_factor
        self.placements = truck_placements(self.deck, placement_step, lane_rule)
        self.fpe = effective_prestress(inp.fpj_ratio, inp.fpu, inp.lump_sum_loss)
        self.girder_ids = self.deck.girder_members.ravel()
        if n_strands is None:
            counts = (Girder1D(inp, n_el, dy_max).design.n_strands,) * self.deck.n_girders
            self.design_rounds = []
            for _ in range(max_design_rounds):
                self._build(counts)
                new = tuple(self._design_count(g) for g in range(self.deck.n_girders))
                self.design_rounds.append(new)
                if new == counts:
                    break
                counts = new
            else:
                raise RuntimeError(f"strand design did not settle: {self.design_rounds}")
        else:
            if len(n_strands) != self.deck.n_girders:
                raise ValueError("n_strands needs one count per girder")
            self.design_rounds = []
            self._build(tuple(int(c) for c in n_strands))
        self._mcr_cache, self._pre_cache, self._w_cache = {}, {}, {}
        self.envelope_el = np.array([self._envelope(g) for g in range(self.deck.n_girders)])

    # ---------------------------------------------------------------- build and design
    def _build(self, counts: tuple) -> None:
        d, inp, mid = self.deck, self.inp, self.deck.mid
        self.n_strands = counts
        self.girders = [Girder1D(replace(inp, girder_spacing=float(w)), self.n_el, self.dy_max, n_strands=c)
                        for w, c in zip(self.widths, counts)]
        self.EI_I = self.EI_fixed.copy()
        for g, gd in enumerate(self.girders):
            self.EI_I[d.girder_members[g]] = gd.EI_I
        self.stiff_I = stiffness(d.model, self.EI_I, self.GJ0)
        sol = solve(self.model_dead, stiff=self.stiff_I)
        self.M_B_dead_mid_I = np.array([line_node_moments(self.model_dead, sol, d.girder_members[g])[mid]
                                        for g in range(d.n_girders)])
        self.M_B_dead_el_I = np.array([midpoint_moments(self.model_dead, sol, d.girder_members[g])
                                       for g in range(d.n_girders)])
        self.M_truck_mid = self._placement_moments()                   # (n_placements, n_girders)
        self.critical = [int(np.argmax(self.M_truck_mid[:, g])) for g in range(d.n_girders)]
        self.M_B_service_mid_I = np.array([self.M_B_dead_mid_I[g] + self.M_truck_mid[self.critical[g], g]
                                           for g in range(d.n_girders)])
        self.soffit_service = [prestress_staged(gd.make_section, gd.is_girder, gd.h_girder, c, inp.strand_area,
                                                self.fpe, float(gd.M_A_nd[mid]),
                                                float(gd.M_A_nd[mid]) + self.M_B_service_mid_I[g]).soffit_stress_service
                               for g, (gd, c) in enumerate(zip(self.girders, counts))]

    def _design_count(self, g: int) -> int:
        gd, mid = self.girders[g], self.deck.mid
        M_A = float(gd.M_A_nd[mid])
        return design_prestress_staged(gd.make_section, gd.is_girder, gd.h_girder, self.inp.strand_area,
                                       self.fpe, M_A, M_A + self.M_B_service_mid_I[g]).n_strands

    def truck_model(self, wheel_lines, axles=None, base=None):
        """Model with the truck set (wheel lines z, axles (x, P)) scaled by (1 + I) x 1.30, plus the
        loads of `base` (default: none). No distribution factor: the grillage distributes."""
        pts, nod = wheel_loads(self.deck, self.axles if axles is None else axles, wheel_lines, self.truck_scale)
        if base is None:
            return with_loads(self.deck.model, q=0.0, point_loads=pts, nodal_loads=nod)
        nodal = dict(base.nodal_loads)
        for dof, v in nod:
            nodal[dof] = nodal.get(dof, 0.0) + v
        return with_loads(base, q=base.udl(), point_loads=list(base.point_loads) + pts,
                          nodal_loads=sorted(nodal.items()))

    def _placement_moments(self) -> np.ndarray:
        """Midspan stage-B moment of every girder for every transverse placement (truck at the
        Stage-1 longitudinal placement), by superposition of one solve per girder line: a wheel
        line at z with lever share s on line g' and arm e adds s/2 (R_g' + e Rt_g')."""
        d, mid = self.deck, self.deck.mid
        R, Rt = [], []
        for g in range(d.n_girders):
            for torque in (False, True):
                model = self._unit_line_model(g, torque)
                sol = solve(model, stiff=self.stiff_I)
                vals = np.array([line_node_moments(model, sol, d.girder_members[h])[mid] for h in range(d.n_girders)])
                (Rt if torque else R).append(vals)
        out = np.zeros((len(self.placements), d.n_girders))
        for p, lines in enumerate(self.placements):
            for zw in lines:
                for g, s, arm in lever_weights(d.z, zw):
                    out[p] += s / 2.0 * (R[g] + arm * Rt[g])
        return self.truck_scale * out

    def _unit_line_model(self, g: int, torque: bool):
        """All axles (full axle load) on girder line g; torque=True: the same loads as torques
        about x with unit arm (consistent nodal split) instead of vertical loads."""
        d = self.deck
        pts, nodal = [], {}
        for x, P in self.axles:
            k, a = d.locate(float(x))
            L = d.x[k + 1] - d.x[k]
            if torque:
                for node, f in ((d.node(g, k), 1.0 - a / L), (d.node(g, k + 1), a / L)):
                    nodal[3 * node + 1] = nodal.get(3 * node + 1, 0.0) + P * f
            else:
                pts.append((int(d.girder_members[g, k]), a, -P))
        return with_loads(d.model, q=0.0, point_loads=pts, nodal_loads=sorted(nodal.items()))

    def _envelope(self, g: int) -> np.ndarray:
        """Total service moment envelope of girder g's elements (midpoints): stage A + stage-B dead
        + the moving truck at g's critical placement (each axle at every node and element midpoint,
        both directions), linear intact grillage."""
        d, inp = self.deck, self.inp
        lines = self.placements[self.critical[g]]
        targets = np.concatenate([d.x, 0.5 * (d.x[:-1] + d.x[1:])])
        best = np.full(self.n_el, -np.inf)
        seen = set()
        for rev in (False, True):
            axles = axle_offsets(inp, rev)
            for off, _ in axles:
                for xt in targets:
                    x0 = round(float(xt - off), 9)
                    if (rev, x0) in seen:
                        continue
                    seen.add((rev, x0))
                    on = [(x0 + o, P) for o, P in axles if -1e-9 <= x0 + o <= inp.span + 1e-9]
                    model = self.truck_model(lines, axles=on)
                    sol = solve(model, stiff=self.stiff_I)
                    best = np.maximum(best, midpoint_moments(model, sol, d.girder_members[g]))
        gd = self.girders[g]
        return gd.M_A_el + self.M_B_dead_el_I[g] + best

    # ---------------------------------------------------------------- nonlinear analysis
    def service_model(self, k: int):
        """Stage-B dead load + the truck at girder k's critical transverse placement."""
        return self.truck_model(self.placements[self.critical[k]], base=self.model_dead)

    def _mcr(self, g: int, eta: float) -> tuple[np.ndarray, np.ndarray]:
        key = (g, float(eta))
        if key not in self._mcr_cache:
            gd = self.girders[g]
            self._mcr_cache[key] = (np.array([gd.mcr(eta, s) for s in gd.sec_el]),
                                    np.array([gd.mcr(eta, s, 0.0) for s in gd.sec_el]))
        return self._mcr_cache[key]

    def precracked(self, g: int, eta: float) -> np.ndarray:
        """Elements of girder g whose service envelope exceeds the sound Mcr(eta)."""
        return self.envelope_el[g] > self._mcr(g, eta)[0]

    def stiffness_fn(self, etas, pre0, sticky: bool = True):
        """EI of every girder element from its stage-B midpoint moment (deck order g, k); the
        element's total moment is M_A + M_B and its reference state the end of stage A (Stage-1
        rule). Elements flagged in `fn.precracked` use the decompression Mcr. sticky=True: an
        element whose total moment exceeds the sound Mcr in any call is flagged from then on.
        sticky=False: the flags are only read and an unflagged element stays uncracked whatever
        its moment; `analyse` flags it once a converged state puts it above the sound Mcr. Both
        modes give the same EI for the same flags and moments at or below the sound Mcr."""
        n, beta = self.n_el, self.inp.beta
        pre = np.array(pre0, dtype=bool).copy()
        mcr = [self._mcr(g, eta) for g, eta in enumerate(etas)]
        e0 = [gd.eps_p0(eta) for gd, eta in zip(self.girders, etas)]

        def fn(M_b: np.ndarray) -> np.ndarray:
            out = np.empty(M_b.size)
            for g, gd in enumerate(self.girders):
                sound, dec = mcr[g]
                for i in range(n):
                    M_ref = float(gd.M_A_el[i])
                    M = M_ref + float(M_b[g * n + i])
                    if sticky and M > sound[i]:
                        pre[g, i] = True
                    if pre[g, i]:
                        m_cr = dec[i]
                    else:                       # sticky: M <= sound here; else uncracked by rule
                        m_cr = sound[i] if sticky else max(sound[i], M)
                    out[g * n + i] = ec2_secant_EI(gd.sec_el[i], e0[g], M, m_cr, beta, M_ref)
            return out

        fn.precracked = pre
        return fn

    def load_path(self, model, base=None, n_steps: int | None = None) -> list:
        """Stage-B loading path ending at `model`: `base` (default the stage-B dead load) first,
        then base + s (model - base) for s = 1/n_steps, ..., 1. One step if model has base's loads."""
        base = self.model_dead if base is None else base
        n_steps = self.n_steps if n_steps is None else int(n_steps)
        if n_steps < 1:
            raise ValueError("n_steps must be >= 1")
        if _same_loads(model, base):
            return [model]
        return [base] + [_blend(base, model, j / n_steps) for j in range(1, n_steps)] + [model]

    def analyse(self, model, etas, gj_factor: float = 1.0, relax: float = 1.0, tol: float = 1e-6,
                tol_EI: float = 1e-6, max_iter: int = 100, n_steps: int | None = None,
                base=None) -> GrillageAnalysis:
        """Secant iteration on the grillage along the load path of `load_path` (stage-B dead load,
        then the rest of `model` in n_steps increments), converged at every step and warm-started
        from the previous one. Inside a run, unflagged elements stay uncracked (no flag is set by
        an intermediate iterate); after each converged run, the unflagged element with the largest
        M/Mcr > 1 (sound Mcr) is flagged pre-cracked and the step is re-converged, one cracking
        event at a time; flags are never cleared.
        max_iter applies to each secant run; a run that does not converge ends the analysis and
        is returned (converged=False). gj_factor scales the GJ of girder elements that are
        cracked (EI below uncracked) - a sensitivity option; 1.0 keeps torsion linear."""
        etas = tuple(float(e) for e in etas)
        ng, n = self.deck.n_girders, self.n_el
        if len(etas) != ng:
            raise ValueError("etas needs one value per girder")
        pre0 = np.array([self.precracked(g, e) for g, e in enumerate(etas)])
        fn = self.stiffness_fn(etas, pre0, sticky=False)
        flags = fn.precracked                                    # read by fn, set here
        sound = np.array([self._mcr(g, e)[0] for g, e in enumerate(etas)])
        M_A = np.array([gd.M_A_el for gd in self.girders])
        ids, EI_I_g = self.girder_ids, self.EI_I[self.girder_ids]

        def runner(step_model):
            def run(EI_g):
                EI = self.EI_fixed.copy()
                EI[ids] = EI_g
                GJ = self.GJ0.copy()
                if gj_factor != 1.0:
                    cracked = EI_g < EI_I_g * (1.0 - 1e-9)
                    GJ[ids[cracked]] *= gj_factor
                sol = solve(step_model, EI, GJ)
                return sol, midpoint_moments(step_model, sol, ids)
            return run

        EI_g = EI_I_g.copy()
        path_max = np.full((ng, n), -np.inf)
        step_its = []
        for step_model in self.load_path(model, base, n_steps):
            its = 0
            while True:
                res = secant_iterate(runner(step_model), fn, EI_g, tol=tol, tol_EI=tol_EI,
                                     max_iter=max_iter, relax=relax)
                its += res.iterations
                if not res.converged:
                    break
                EI_g = res.EI
                total = M_A + res.M.reshape(ng, n)
                ratio = np.where(flags, -np.inf, total / sound)
                if ratio.max() <= 1.0:                           # accepted: no unflagged element above Mcr
                    path_max = np.maximum(path_max, total)
                    break
                # one cracking event at a time: the element furthest above its sound Mcr cracks
                # first; the others are re-checked after the load it sheds has been redistributed
                ev = np.unravel_index(int(np.argmax(ratio)), ratio.shape)
                path_max[ev] = max(path_max[ev], total[ev])
                flags[ev] = True                                 # fn sees the new flag
            step_its.append(its)
            if not res.converged:
                break
        d, mid = self.deck, self.deck.mid
        sol = res.solution
        M_mid = np.array([line_node_moments(step_model, sol, d.girder_members[g])[mid] for g in range(ng)])
        defl = np.array([-sol.D[3 * d.node(g, mid)] for g in range(ng)])
        return GrillageAnalysis(res, etas, flags.copy(), M_mid, defl, tuple(step_its), path_max)

    def etas(self, k: int, eta: float) -> tuple:
        e = [1.0] * self.deck.n_girders
        e[k] = float(eta)
        return tuple(e)

    def w_pred_detail(self, eta: float, k: int, gj_factor: float = 1.0) -> dict:
        """Crack width at Gk's midspan node under the total dead load (the inspection state):
        stage A from statics, stage B through the nonlinear grillage."""
        key = (float(eta), int(k), float(gj_factor))
        if key not in self._w_cache:
            an = self.analyse(self.model_dead, self.etas(k, eta), gj_factor)
            gd, mid = self.girders[k], self.deck.mid
            M = float(gd.M_A_nd[mid]) + float(an.M_B_mid[k])
            w = gd.crack_at(eta, M, gd.sec_nd[mid])["width"] if an.result.converged else float("nan")
            self._w_cache[key] = {"width": w, "M_total": M, "converged": an.result.converged,
                                  "iterations": an.result.iterations}
        return self._w_cache[key]

    def w_pred(self, eta: float, k: int, gj_factor: float = 1.0) -> float:
        """Width in metres; NaN if the secant loop did not converge (the inverse refuses NaN)."""
        return self.w_pred_detail(eta, k, gj_factor)["width"]

    def w_pred_fn(self, k: int):
        return lambda eta: self.w_pred(eta, k)

    def indicator(self, eta: float, k: int, gj_factor: float = 1.0, relax: float = 1.0) -> dict:
        d, mid = self.deck, self.deck.mid
        serv_model = self.service_model(k)
        intact = {"dead": self.analyse(self.model_dead, self.etas(k, 1.0)),
                  "service": self.analyse(serv_model, self.etas(k, 1.0))}
        dmg = {"dead": self.analyse(self.model_dead, self.etas(k, eta), gj_factor, relax),
               "service": self.analyse(serv_model, self.etas(k, eta), gj_factor, relax)}
        Mi, Md = intact["service"].M_B_mid, dmg["service"].M_B_mid
        truck_i = Mi - intact["dead"].M_B_mid
        truck_d = Md - dmg["dead"].M_B_mid
        s = dmg["service"].result
        EI_g = s.EI.reshape(d.n_girders, self.n_el)
        EI_I_g = self.EI_I[self.girder_ids].reshape(d.n_girders, self.n_el)
        return {
            "eta": float(eta), "k": int(k), "gj_factor": gj_factor,
            "w_pred_mm": self.w_pred(eta, k, gj_factor) * 1e3,
            "M_B_service_mid_kNm": {"intact": (Mi / 1e3).tolist(), "damaged": (Md / 1e3).tolist()},
            "M_truck_mid_kNm": {"intact": (truck_i / 1e3).tolist(), "damaged": (truck_d / 1e3).tolist()},
            "share_service": {"intact": (Mi / Mi.sum()).tolist(), "damaged": (Md / Md.sum()).tolist()},
            "shed_percent_of_Gk": ((Md - Mi) / Mi[k] * 100.0).tolist(),
            "deflection_service_mm": {"intact": (intact["service"].deflection_mid * 1e3).tolist(),
                                      "damaged": (dmg["service"].deflection_mid * 1e3).tolist()},
            "deflection_truck_mm": {
                "intact": ((intact["service"].deflection_mid - intact["dead"].deflection_mid) * 1e3).tolist(),
                "damaged": ((dmg["service"].deflection_mid - dmg["dead"].deflection_mid) * 1e3).tolist()},
            "deflection_limit_mm": self.inp.span / self.inp.deflection_limit_ratio * 1e3,
            "EI_ratio_Gk": (EI_g[k] / EI_I_g[k]).tolist(),
            "cracked_elements": [int(np.sum(EI_g[g] < EI_I_g[g] * (1.0 - 1e-9))) for g in range(d.n_girders)],
            "precracked_Gk": dmg["service"].precracked[k].tolist(),
            "x_mid": (0.5 * (d.x[:-1] + d.x[1:])).tolist(),
            "converged": {"dead": dmg["dead"].result.converged, "service": s.converged},
            # secant updates of the last load step; passes_per_step lists every step (dead load first)
            "iterations": {"dead": dmg["dead"].result.iterations, "service": s.iterations},
            "passes_per_step": {"dead": list(dmg["dead"].step_iterations),
                                "service": list(dmg["service"].step_iterations)},
            "precracked_count": dmg["service"].precracked.sum(axis=1).tolist(),
            "history_service": list(s.history), "ei_residual_service": list(s.ei_residual),
            "history_dead": list(dmg["dead"].result.history),
        }
