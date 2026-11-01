"""Stage 1: the whole chain on one simply supported girder (Gk), built in two stages.

Stage A: the bare girder carries the prestress and the weight of the girder and the wet slab
(moments from statics). Stage B: wearing surface, barriers and the truck act on the composite
girder through the FE model. Each element and node has its own staged section, because the
slab's locked-in strain depends on the local stage-A moment. eta scales today's strand strain;
the slab strain locked in at construction (eta = 1) does not change.

Crack width at rest (total dead load) -> eta through the inverse; stiffness per element from
its own total moment, measured from the end of stage A; deflection under stage-B dead load +
truck. Elements whose total service moment under the moving-truck envelope has exceeded Mcr(eta)
are 'pre-cracked' (traffic cracked them): their tensile strength is taken as zero, so their effective Mcr is the
decompression moment.
"""
import numpy as np

from .crackwidth import EC2Crack, effective_tension_height
from .geometry import area_below, discretize, girder_outline, gross_properties, slab, top_of
from .loads import (impact_factor, max_truck_moment_at, midspan_moment_statics, service_truck_loads,
                    stage_a_line, stage_b_line, truck_fraction)
from .materials import ec2_ecm, ec2_fctm
from .model import simply_supported, with_loads
from .nonlinear import solve_secant
from .params import Inputs
from .prestress import effective_prestress
from .section import PCSection, cracking_moment, delta_sigma_p, solve_state, uncracked_EI
from .solve import midpoint_moments, node_moments, solve
from .staging import design_prestress_staged, girder_only, prestress_staged, staged_section
from .stiffness import branson_secant_EI, ec2_kappa, ec2_secant_EI


class Girder1D:
    def __init__(self, inp: Inputs = Inputs(), n_el: int = 40, dy_max: float = 1e-3,
                 n_strands: int | None = None):
        """n_strands=None designs the strand count (smallest with no soffit tension under service);
        an integer analyses that as-built count instead (e.g. an existing girder whose inputs are
        perturbed in a sensitivity study does not get new strands)."""
        if n_el % 2:
            raise ValueError("n_el must be even so that midspan is a node")
        self.inp, self.n_el = inp, n_el
        L = inp.span
        self.Ec_girder = ec2_ecm(inp.fck_girder)
        self.fct_girder = ec2_fctm(inp.fck_girder)
        self.girder_parts = girder_outline(inp.girder_outline_in, self.Ec_girder, self.fct_girder)
        self.h_girder = top_of(self.girder_parts)
        self.parts = self.girder_parts + [slab(self.h_girder, inp.slab_thickness, inp.girder_spacing,
                                               ec2_ecm(inp.fck_slab), ec2_fctm(inp.fck_slab))]
        self.h = top_of(self.parts)
        self.fibers = discretize(self.parts, dy_max)
        self.is_girder = self.fibers.y < self.h_girder
        self.q_a = stage_a_line(inp, gross_properties(self.girder_parts)[0])
        self.q_b = stage_b_line(inp)
        self.truck = service_truck_loads(inp)
        base = simply_supported(L, n_el)
        self.mid_node = n_el // 2

        def stage_a(x):                     # statics: the bare girder is simply supported
            return self.q_a * x * (L - x) / 2.0

        self.M_A_el, self.M_A_nd = stage_a(base.midpoints), stage_a(base.x)
        self.M_service = midspan_moment_statics(L, self.q_a + self.q_b, self.truck)
        fpe = effective_prestress(inp.fpj_ratio, inp.fpu, inp.lump_sum_loss)
        M_A_mid = float(self.M_A_nd[self.mid_node])
        if n_strands is None:
            self.design = design_prestress_staged(self.make_section, self.is_girder, self.h_girder,
                                                  inp.strand_area, fpe, M_A_mid, self.M_service)
        else:
            self.design = prestress_staged(self.make_section, self.is_girder, self.h_girder, n_strands,
                                           inp.strand_area, fpe, M_A_mid, self.M_service)
        self.composite = self.make_section(self.design.A_p)
        e_c = self.design.eps_p0            # construction state: eta = 1
        self.sec_el = [staged_section(self.composite, self.is_girder, self.h_girder, e_c, float(m))
                       for m in self.M_A_el]
        self.sec_nd = [staged_section(self.composite, self.is_girder, self.h_girder, e_c, float(m))
                       for m in self.M_A_nd]
        self.EI_I = uncracked_EI(self.composite)
        self.EI_girder = uncracked_EI(girder_only(self.composite, self.is_girder, self.h_girder))
        self.crack = EC2Crack(inp.k1, inp.k2, inp.k3, inp.k4, inp.kt, inp.bond_xi, inp.Ep,
                              self.Ec_girder, self.fct_girder, inp.strand_cover,
                              inp.strand_diameter, self.design.A_p)
        self.model_dead = with_loads(base, q=-self.q_b)
        self.model_service = with_loads(base, q=-self.q_b,
                                        point_loads=[(x, -P) for x, P in self.truck])
        # Determinate girder: stage-B moments do not depend on stiffness (test_nonlinear proves the
        # secant loop stops after one pass), so the total moments are computed once.
        sol_dead = solve(self.model_dead, self.EI_I)
        self.M_dead_nd = self.M_A_nd + node_moments(self.model_dead, sol_dead)
        self.M_dead_mid = float(self.M_dead_nd[self.mid_node])
        self.M_service_el = self.M_A_el + midpoint_moments(self.model_service,
                                                           solve(self.model_service, self.EI_I))
        # Moving-truck envelope (each element's own worst truck position): decides which elements
        # traffic has cracked. Deflections still use the midspan placement (self.truck).
        scale = truck_fraction(inp) * (1.0 + impact_factor(inp)) * inp.live_load_factor
        self.M_service_env_el = (self.M_A_el + midpoint_moments(self.model_dead, sol_dead)
                                 + scale * np.array([max_truck_moment_at(inp, float(x))
                                                     for x in base.midpoints]))

    def make_section(self, A_p: float) -> PCSection:
        return PCSection(self.fibers, self.inp.strand_centroid, A_p, self.inp.Ep, self.h)

    def eps_p0(self, eta: float) -> float:
        return eta * self.design.eps_p0

    def mcr(self, eta: float, sec: PCSection, fct: float | None = None) -> float:
        """Total moment at which this staged section's soffit reaches fct (0 = decompression)."""
        return cracking_moment(sec, self.eps_p0(eta), self.fct_girder if fct is None else fct)

    def crack_at(self, eta: float, M: float, sec: PCSection) -> dict:
        st = solve_state(sec, self.eps_p0(eta), M, cracked=True)
        ds = delta_sigma_p(sec, st)
        if ds <= 0.0:
            return {"delta_sigma_p": ds, "Ac_eff": None, "width": 0.0, "neutral_axis": None}
        y_na = st.neutral_axis if st.kappa != 0.0 else float("nan")
        if not (0.0 < y_na < self.h_girder):   # compression zone ends in the slab
            y_na = float(sec.fib.y[(st.strain(sec.fib.y) - sec.free_strain()) < 0.0].min())
        hce = effective_tension_height(self.h, sec.y_p, y_na)
        Ac = area_below(self.girder_parts, hce)
        return {"delta_sigma_p": ds, "Ac_eff": Ac, "width": self.crack.width(ds, Ac),
                "neutral_axis": y_na}

    def w_pred(self, eta: float) -> float:
        """Predicted crack width (m) at midspan under the total dead load (the inspection state)."""
        return self.crack_at(eta, self.M_dead_mid, self.sec_nd[self.mid_node])["width"]

    def precracked_zone(self, eta: float) -> np.ndarray:
        """Elements whose total service moment under the moving-truck envelope has exceeded Mcr(eta)."""
        return np.array([m > self.mcr(eta, s) for m, s in zip(self.M_service_env_el, self.sec_el)])

    def damaged_zone(self, length: float) -> np.ndarray:
        """Elements whose midpoints lie within length/2 of midspan (exact when length is a
        multiple of two element lengths)."""
        return np.abs(self.model_service.midpoints - self.inp.span / 2.0) < length / 2.0

    def stiffness_fn(self, eta: float, precracked: np.ndarray, method: str = "ec2",
                     damaged: np.ndarray | None = None):
        """EI per element from its stage-B moment (what the FE model returns); the element's
        total moment is M_A + M_B and its reference state is the end of stage A.
        damaged: elements allowed to crack (None = all); the others keep their uncracked EI, as in
        the proposal's damage confined to one element."""
        e0 = self.eps_p0(eta)
        mcr = [self.mcr(eta, s, 0.0 if pre else None) for s, pre in zip(self.sec_el, precracked)]

        def fn(M_b: np.ndarray) -> np.ndarray:
            out = np.empty(M_b.size)
            for i, mb in enumerate(M_b):
                if damaged is not None and not damaged[i]:
                    out[i] = uncracked_EI(self.sec_el[i])
                    continue
                M_ref = float(self.M_A_el[i])
                M = M_ref + float(mb)
                if method == "ec2":
                    out[i] = ec2_secant_EI(self.sec_el[i], e0, M, mcr[i], self.inp.beta, M_ref)
                else:
                    out[i] = branson_secant_EI(self.sec_el[i], e0, M, mcr[i], M_ref)
            return out

        return fn

    def analyse(self, eta: float, model, method: str = "ec2", damaged: np.ndarray | None = None):
        return solve_secant(model, self.stiffness_fn(eta, self.precracked_zone(eta), method, damaged), self.EI_I)

    def kappa_curves(self, eta: float, M_grid) -> dict:
        """Midspan moment-curvature for the report figure: total moment on the x axis; curvature
        measured from the end of stage A (M = M_A), as the FE's secant stiffness is. Each EC2 curve
        is measured from its own ec2_kappa(M_A) (the FE's reference); the uncracked and fully
        cracked curves from the uncracked kappa(M_A). The fully cracked curve is drawn only above
        the decompression moment (below it, prestress alone would crack the top)."""
        sec = self.sec_nd[self.mid_node]
        e0 = self.eps_p0(eta)
        beta = self.inp.beta
        M_A = float(self.M_A_nd[self.mid_node])
        mcr_sound, mcr_pre = self.mcr(eta, sec), self.mcr(eta, sec, 0.0)
        k_ref = solve_state(sec, e0, M_A, cracked=False).kappa
        k_ref_sound = ec2_kappa(sec, e0, M_A, mcr_sound, beta)
        k_ref_pre = ec2_kappa(sec, e0, M_A, mcr_pre, beta)
        out = {"uncracked": [], "cracked": [], "ec2_sound": [], "ec2_precracked": []}
        for M in M_grid:
            M = float(M)
            out["uncracked"].append(solve_state(sec, e0, M, cracked=False).kappa - k_ref)
            out["cracked"].append(solve_state(sec, e0, M, cracked=True).kappa - k_ref
                                  if M >= mcr_pre else float("nan"))
            out["ec2_sound"].append(ec2_kappa(sec, e0, M, mcr_sound, beta) - k_ref_sound)
            out["ec2_precracked"].append(ec2_kappa(sec, e0, M, mcr_pre, beta) - k_ref_pre)
        return out

    def visible_crack_zone(self, eta: float) -> float:
        """Length of span (m) whose nodes show an open crack at rest. On synthetic data this is a
        model output, not a held-out check: the masks come from the same model."""
        open_ = [self.crack_at(eta, float(M), s)["width"] > 0.0 for M, s in zip(self.M_dead_nd, self.sec_nd)]
        xs = self.model_dead.x[open_]
        return float(xs.max() - xs.min()) if xs.size else 0.0

    def indicator(self, eta: float, method: str = "ec2", damaged: np.ndarray | None = None) -> dict:
        dead = self.analyse(eta, self.model_dead, method, damaged)
        serv = self.analyse(eta, self.model_service, method, damaged)
        v_dead = -dead.solution.D[2 * self.mid_node]
        v_serv = -serv.solution.D[2 * self.mid_node]
        L = self.inp.span
        sec_mid = self.sec_nd[self.mid_node]
        k_transfer = solve_state(girder_only(self.composite, self.is_girder, self.h_girder),
                                 self.design.eps_p0, 0.0, cracked=False).kappa
        return {
            "eta": eta,
            "method": method,
            "w_pred_mm": self.w_pred(eta) * 1e3,
            "Mcr_mid_kNm": self.mcr(eta, sec_mid) / 1e3,
            "Mdec_mid_kNm": self.mcr(eta, sec_mid, 0.0) / 1e3,
            "M_dead_mid_kNm": self.M_dead_mid / 1e3,
            "camber_mm": -k_transfer * L * L / 8.0 * 1e3,            # prestress-only camber of the bare girder at fpe (eta = 1), no self-weight (independent of eta)
            "deflection_stage_a_mm": 5.0 * self.q_a * L**4 / (384.0 * self.EI_girder) * 1e3,
            "deflection_dead_mm": v_dead * 1e3,                      # stage-B dead load, composite
            "deflection_service_mm": v_serv * 1e3,                   # stage-B dead load + truck
            "deflection_truck_mm": (v_serv - v_dead) * 1e3,
            "deflection_limit_mm": L / self.inp.deflection_limit_ratio * 1e3,
            "visible_crack_zone_m": self.visible_crack_zone(eta),
            "EI_ratio": (serv.EI / self.EI_I).tolist(),
            "x_mid": self.model_service.midpoints.tolist(),
            "x_nodes": self.model_service.x.tolist(),
            "v_service_mm": (-serv.solution.D[0::2] * 1e3).tolist(),
            "precracked": self.precracked_zone(eta).tolist(),
            "iterations_service": serv.iterations,
        }
