"""Regenerate every Stage-1 result and figure from committed inputs.

Reads data/claimbridge/*_spike.json (thesis-pipeline output) and writes
results/stage1.json + results/figures/*.png. Does not need the thesis pipeline.
"""
import json
import math
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem import figures as F  # noqa: E402
from crackfem.claimbridge import governing_width, refused_below_resolution  # noqa: E402
from crackfem.girder1d import Girder1D  # noqa: E402
from crackfem.inverse import estimate_eta, eta_decompression, is_monotone_nonincreasing  # noqa: E402
from crackfem.params import Inputs  # noqa: E402
from crackfem.solve import solve  # noqa: E402
from crackfem.uq import DEFAULT_SEED, monte_carlo_eta, tornado  # noqa: E402

MC_N = 2000
MESH_NS = (10, 20, 40, 80, 160)
TORNADO = {"girder_spacing": (1.8, 2.6), "fck_girder": (35e6, 45e6), "strand_cover": (0.035, 0.055),
           "strand_centroid": (0.08, 0.12), "slab_thickness": (0.18, 0.22), "bond_xi": (0.5, 0.7)}
TORNADO_LABELS = {"girder_spacing": ("girder spacing", 1.0, "m"), "fck_girder": ("girder fck", 1e-6, "MPa"),
                  "strand_cover": ("strand cover", 1e3, "mm"), "strand_centroid": ("strand centroid", 1e3, "mm"),
                  "slab_thickness": ("slab thickness", 1e3, "mm"), "bond_xi": ("bond factor ξ", 1.0, "")}
GATE_PX = 3.0          # ClaimBridge refuses widths below 3 px (min measurable width = 3 GSD)
SIGMA_PX = 1.0         # measurement sigma = u/2 with u = 2 GSD (the Monte Carlo convention, uq.py)
Z95_ONE_SIDED = 1.645
# Stage-2 (grillage) input fields, Plan 2 Task 1: not part of the Stage-1 results, so a rerun
# of this script leaves stage1.json["inputs"] as it was before Stage 2.
STAGE2_INPUTS = ("overhang_width_ft", "barrier_width_ft", "wheel_spacing_ft", "truck_clearance_ft",
                 "wheel_curb_clearance_ft", "lane_width_ft", "poisson_concrete", "diaphragm_depth",
                 "diaphragm_thickness")


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return jsonable(o.item())
    if isinstance(o, float) and not math.isfinite(o):
        return None                   # strict JSON has no NaN/Infinity
    return o


def run_case(g: Girder1D, rec: dict) -> dict:
    c = rec["case"]
    gov = governing_width(rec["parsed"])
    out = {"name": c["name"], "gsd_mm": c["gsd_mm"], "eta_true": c["eta_true"], "w_true_mm": c["w_true_mm"],
           "width_px": c["width_px"]}
    if gov is None:
        gates = sorted({k["gate"] for k in rec["parsed"]["cracks"] if k["gate"] != "ok"})
        if not refused_below_resolution(rec["parsed"]):
            # An unresolved crack ("no-valid-points", an unknown gate) or a report without cracks: the
            # governing width may exceed the gate, so there is no bound on eta (field-data safety).
            return {**out, "status": "unresolved", "eta_point": None, "eta_mc": None, "reason": gates,
                    "eta_lower_bound_from_gate": None, "gate_inverse_status": None,
                    "eta_lower_bound_95": None, "gate95_inverse_status": None}
        # Every crack is below the gate, so the measured width is below GATE_PX and eta lies above the
        # eta whose predicted width is GATE_PX px (nominal bound). With measurement sigma = 1 px the true
        # width is below GATE_PX + 1.645 px at one-sided 95 %, which gives the lower (95 %) bound.
        gsd_m = c["gsd_mm"] * 1e-3
        gate = estimate_eta(g.w_pred, GATE_PX * gsd_m)
        gate95 = estimate_eta(g.w_pred, (GATE_PX + Z95_ONE_SIDED * SIGMA_PX) * gsd_m)
        return {**out, "status": "refused", "eta_point": None, "eta_mc": None,
                "reason": gates,
                "eta_lower_bound_from_gate": gate.eta if gate.status == "ok" else None,
                "gate_inverse_status": gate.status,
                "eta_lower_bound_95": gate95.eta if gate95.status == "ok" else None,
                "gate95_inverse_status": gate95.status}
    w, u = gov["max_width_mm"] * 1e-3, gov["uncertainty_mm"] * 1e-3
    est = estimate_eta(g.w_pred, w)
    t0 = time.perf_counter()
    mc = monte_carlo_eta(g.w_pred, w, u, n=MC_N)
    t_mc = time.perf_counter() - t0
    s = mc.summary()
    band = {k: g.indicator(s[k]) for k in ("p05", "p50", "p95")} if s else {}
    # Deflection falls as eta rises, so its 5th percentile comes from eta's 95th and vice versa.
    from_eta = {"p05": "p95", "p50": "p50", "p95": "p05"}

    def dband(key):
        return {**{k: band[v][key] for k, v in from_eta.items()}, "from_eta_percentile": from_eta} if band else None

    return {**out, "status": "measured", "w_meas_mm": gov["max_width_mm"], "u95_mm": gov["uncertainty_mm"],
            "rating": gov["rating"], "rating_label": gov["rating_label"],
            "eta_point": est.eta, "eta_status": est.status, "eta_mc": s, "mc_status_counts": mc.status_counts,
            "mc_seconds": t_mc, "mc_samples": mc.etas,
            "deflection_service_mm_at": {k: v["deflection_service_mm"] for k, v in band.items()},
            "deflection_truck_mm_at": {k: v["deflection_truck_mm"] for k, v in band.items()},
            "deflection_service_band_mm": dband("deflection_service_mm"),
            "deflection_truck_band_mm": dband("deflection_truck_mm")}


def tornado_label(name: str, lo: float, hi: float) -> str:
    text, scale, unit = TORNADO_LABELS.get(name, (name, 1.0, ""))
    return f"{text} {lo * scale:g}–{hi * scale:g} {unit}".strip()


def stage1_inputs() -> dict:
    """The Inputs fields the Stage-1 results depend on (Stage-2 fields left out)."""
    return {k: v for k, v in asdict(Inputs()).items() if k not in STAGE2_INPUTS}


def main() -> None:
    res: dict = {}
    g = Girder1D()
    res["inputs"] = stage1_inputs()
    res["design"] = asdict(g.design)
    res["girder"] = {"q_stage_a_kN_m": g.q_a / 1e3, "q_stage_b_kN_m": g.q_b / 1e3,
                     "M_A_mid_kNm": float(g.M_A_nd[g.mid_node]) / 1e3, "M_dead_mid_kNm": g.M_dead_mid / 1e3,
                     "M_service_kNm": g.M_service / 1e3, "EI_I_GNm2": g.EI_I / 1e9,
                     "EI_girder_GNm2": g.EI_girder / 1e9, "truck": g.truck, "h_m": g.h}
    grid = np.linspace(0.02, 1.0, 99)
    res["w_curve"] = {"eta": grid, "w_mm": [g.w_pred(float(e)) * 1e3 for e in grid],
                      "monotone": is_monotone_nonincreasing(g.w_pred, grid)}
    res["eta_decompression"] = eta_decompression(g.w_pred)

    data = ROOT / "data" / "claimbridge"
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(data.glob("*_spike.json"))]
    res["cases"] = [run_case(g, r) for r in recs]

    measured = [c for c in res["cases"] if c["status"] == "measured"]
    best = min(measured, key=lambda c: (c["gsd_mm"], -c["eta_true"]))
    eta_ref = best["eta_true"]            # finest GSD, least damage: the first detectable case
    res["eta_ref"] = eta_ref
    res["indicator_intact"] = g.indicator(1.0)
    res["indicator_ref"] = g.indicator(eta_ref)
    res["indicator_ref_branson"] = g.indicator(eta_ref, method="branson")
    M_A_mid = float(g.M_A_nd[g.mid_node])
    Mcr_mid = g.mcr(eta_ref, g.sec_nd[g.mid_node])               # sound section, eta_ref; the exact float
    # kappa_curves compares against (a kN·m round trip can land one ulp on the cracked side)
    M_grid = np.linspace(M_A_mid, 1.3 * g.M_service, 60)
    if M_A_mid < Mcr_mid < M_grid[-1]:                            # resolve the EC2 zeta jump at Mcr exactly
        M_grid = np.unique(np.concatenate([M_grid, [Mcr_mid, Mcr_mid * (1.0 + 1e-9)]]))
    res["mkappa"] = {"M_kNm": M_grid / 1e3, **g.kappa_curves(eta_ref, M_grid)}

    res["mesh"] = {"N": list(MESH_NS),
                   "intact_mm": [Girder1D(n_el=n).indicator(1.0)["deflection_service_mm"] for n in MESH_NS],
                   "damaged_mm": [Girder1D(n_el=n).indicator(eta_ref)["deflection_service_mm"] for n in MESH_NS]}

    w_best, u_best = best["w_meas_mm"] * 1e-3, best["u95_mm"] * 1e-3
    res["mc_stability"] = {str(n): monte_carlo_eta(g.w_pred, w_best, u_best, n=n).summary() for n in (1000, 10000)}
    res["mc_stability"]["2000_seed_plus_1"] = monte_carlo_eta(g.w_pred, w_best, u_best, n=MC_N,
                                                              seed=DEFAULT_SEED + 1).summary()

    def eta_of(girder: Girder1D) -> float:
        est = estimate_eta(girder.w_pred, w_best)
        return est.eta if est.status == "ok" else float("nan")

    # The bridge exists: a perturbed input does not give the girder new strands.
    n_built = g.design.n_strands
    res["tornado"] = {"case": best["name"], "strand_count": f"as built ({n_built}, held)",
                      "rows": tornado(lambda inp: eta_of(Girder1D(inp, n_strands=n_built)), Inputs(), TORNADO)}
    res["tornado_redesigned"] = {"case": best["name"], "strand_count": "re-designed for each input (record only)",
                                 "rows": tornado(lambda inp: eta_of(Girder1D(inp)), Inputs(), TORNADO)}

    t0 = time.perf_counter(); solve(g.model_service, g.EI_I); t_lin = time.perf_counter() - t0
    t0 = time.perf_counter(); g.analyse(eta_ref, g.model_service); t_nl = time.perf_counter() - t0
    t0 = time.perf_counter(); estimate_eta(g.w_pred, w_best); t_inv = time.perf_counter() - t0
    res["timing_s"] = {"linear_solve": t_lin, "nonlinear_analysis": t_nl, "inverse": t_inv,
                       "monte_carlo_best_case": best["mc_seconds"], "mc_samples": MC_N}

    out = ROOT / "results"
    figs = out / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    (out / "stage1.json").write_text(json.dumps(jsonable(res), indent=2, ensure_ascii=False, allow_nan=False),
                                     encoding="utf-8")

    w_max = max(res["w_curve"]["w_mm"])
    gsds = sorted({c["gsd_mm"] for c in res["cases"] if GATE_PX * c["gsd_mm"] <= w_max})
    F.fig_w_eta(res["w_curve"]["eta"], res["w_curve"]["w_mm"], measured, figs / "fig1_w_eta.png",
                eta_dec=res["eta_decompression"],
                gate_widths_mm={f"{GATE_PX:g} px at GSD {s} mm/px = {GATE_PX * s:.2f} mm": GATE_PX * s for s in gsds})
    mk, sound = res["mkappa"], "EC2 7.4.3 – sound (fct = fctm)"
    F.fig_mkappa(mk["M_kNm"], {"uncracked": mk["uncracked"],
                               "fully cracked (below decompression: not shown)": mk["cracked"],
                               sound: mk["ec2_sound"],
                               "EC2 7.4.3 – pre-cracked (fct = 0, used by the FE here)": mk["ec2_precracked"]},
                 figs / "fig2_mkappa.png", title=f"Midspan section, η = {eta_ref}",
                 break_at={sound: Mcr_mid / 1e3},
                 hlines={f"M dead, at rest = {g.M_dead_mid / 1e3:.0f} kN·m": g.M_dead_mid / 1e3,
                         f"M service, dead + truck = {g.M_service / 1e3:.0f} kN·m": g.M_service / 1e3})
    F.fig_ei_ratio(res["indicator_ref"]["x_mid"], {"η = 1": res["indicator_intact"]["EI_ratio"],
                                                   f"η = {eta_ref}": res["indicator_ref"]["EI_ratio"]},
                   figs / "fig3_ei_ratio.png")
    ii, ir = res["indicator_intact"], res["indicator_ref"]
    F.fig_deflection(ir["x_nodes"], {"η = 1": ii["v_service_mm"], f"η = {eta_ref}": ir["v_service_mm"]},
                     figs / "fig4_deflection.png",
                     title=(f"Truck only, midspan: {ii['deflection_truck_mm']:.1f} mm (η = 1), "
                            f"{ir['deflection_truck_mm']:.1f} mm (η = {eta_ref}); "
                            f"limit L/{g.inp.deflection_limit_ratio:g} = {ir['deflection_limit_mm']:.1f} mm"))
    F.fig_eta_hist(best["mc_samples"], best["eta_true"], figs / "fig5_eta_hist.png",
                   title=(f"{best['name']}: GSD {best['gsd_mm']} mm/px, w = {best['w_meas_mm']} ± "
                          f"{best['u95_mm']} mm; measurement noise only"),
                   p05=best["eta_mc"]["p05"], p95=best["eta_mc"]["p95"])
    intact = res["mesh"]["intact_mm"]
    spread = (max(intact) - min(intact)) / abs(intact[-1])
    F.fig_mesh(res["mesh"]["N"], {f"η = {eta_ref}": res["mesh"]["damaged_mm"]}, figs / "fig6_mesh.png",
               note=f"intact (η = 1, uniform EI): exact at the nodes for every N (spread {spread:.0e}, round-off)")
    F.fig_cases(res["cases"], figs / "fig7_cases.png")
    F.fig_tornado(res["tornado"]["rows"], figs / "fig8_tornado.png",
                  labels={r["name"]: tornado_label(r["name"], r["low"], r["high"]) for r in res["tornado"]["rows"]},
                  title=f"η sensitivity, as-built {n_built}-strand girder, case {best['name']}")

    print(f"eta_decompression = {res['eta_decompression']:.4f}")
    for c in res["cases"]:
        print(c["name"], c["status"], c.get("eta_point"), c.get("eta_mc"))
    for key in ("tornado", "tornado_redesigned"):
        print(key, {r["name"]: round(r["swing"], 4) for r in res[key]["rows"]})
    print("timing", res["timing_s"])


if __name__ == "__main__":
    main()
