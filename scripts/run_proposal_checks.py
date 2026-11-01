"""Close the proposal's Expected Outcomes 2 and 3 and record one beam-element hand calculation.

Writes results/proposal_checks.json and figures 15-17 with the models of run_all.py (Stage 1)
and run_stage2.py (grillage). It does not change stage1.json or stage2.json.

- Expected Outcome 2: truck deflection and stiffness ratio against the measured crack width.
- Expected Outcome 3: damage confined to the elements at midspan (the proposal's one damaged
  element, length tied to the mesh) against a fixed physical length and the distributed model.
- Hand calculation: element stiffness, consistent loads, DOF map and equilibrium, N = 10.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem import figures as F  # noqa: E402
from crackfem.elements import beam_k, udl_equivalent  # noqa: E402
from crackfem.girder1d import Girder1D  # noqa: E402
from crackfem.grillage import Grillage  # noqa: E402
from crackfem.inverse import estimate_eta, eta_decompression  # noqa: E402
from crackfem.solve import solve  # noqa: E402

ETA_REF = 0.36                      # Stage-1 reference case (case2), as in run_all.py
MESH_NS = (10, 20, 40, 80, 160)
LENGTH_STEPS = (1, 2, 4, 8, 16, 24, 32, 48, 64, 80)   # damaged length = steps x 2 element lengths, N = 160
N_EL_GRID = 20                      # grillage mesh of run_stage2.py
GIRDERS = {"G1 exterior": 0, "G3 interior": 2}
GATE_PX = 3.0
GSDS_MM = (0.02, 0.05, 0.10)
BAND_GSD_MM = 0.02                  # measured cases drawn with their band (finest GSD)
N_ETA_STAGE1, N_ETA_GRID = 25, 13
ETA_MIN = 0.05


def midspan_ratio(ratios) -> float:
    """EI ratio at midspan: the smaller of the two elements that meet at the midspan node."""
    n = len(ratios)
    return float(min(ratios[n // 2 - 1], ratios[n // 2]))


def stage1_rows(g: Girder1D, etas) -> list:
    """alpha_equivalent = truck deflection at eta = 1 over truck deflection at eta: the factor a
    uniform EI would need to give the same truck deflection (the proposal's Ie / Ig, globally)."""
    base = g.indicator(1.0)["deflection_truck_mm"]
    rows = []
    for eta in etas:
        ind = g.indicator(float(eta))
        rows.append({"eta": float(eta), "w_mm": ind["w_pred_mm"],
                     "deflection_truck_mm": ind["deflection_truck_mm"],
                     "EI_ratio_midspan": midspan_ratio(ind["EI_ratio"]), "EI_ratio_min": min(ind["EI_ratio"]),
                     "alpha_equivalent": base / ind["deflection_truck_mm"]})
    return rows


def stage1_case_band(g: Girder1D, case: dict, base: float) -> dict:
    """A measured case carried through: the response at the Monte Carlo 5/50/95 % eta, sorted
    low to high (deflection rises as eta falls, the stiffness ratios fall)."""
    ind = {q: g.indicator(case["eta_mc"][q]) for q in ("p05", "p50", "p95")}
    out = {"name": case["name"], "w": case["w_meas_mm"], "u": case["u95_mm"]}
    for key, value in (("deflection_truck_mm", lambda i: i["deflection_truck_mm"]),
                       ("alpha_equivalent", lambda i: base / i["deflection_truck_mm"]),
                       ("EI_ratio_midspan", lambda i: midspan_ratio(i["EI_ratio"]))):
        v = sorted(value(ind[q]) for q in ("p05", "p95"))
        out[key] = {"lo": v[0], "mid": value(ind["p50"]), "hi": v[1]}
    return out


def end_point_rows(g: Girder1D, cases) -> list:
    """The proposal's band: eta solved again at w + u and w - u (u = 2·GSD), beside the
    Monte Carlo 5-95 % interval (sigma = u/2) that the report uses instead."""
    rows = []
    for c in cases:
        w, u = c["w_meas_mm"] * 1e-3, c["u95_mm"] * 1e-3
        wide, narrow = estimate_eta(g.w_pred, w + u), estimate_eta(g.w_pred, w - u)   # wider crack -> lower eta
        rows.append({"name": c["name"], "gsd_mm": c["gsd_mm"], "w_mm": c["w_meas_mm"], "u_mm": c["u95_mm"],
                     "eta_point": c["eta_point"], "eta_at_w_plus_u": wide.eta, "eta_at_w_minus_u": narrow.eta,
                     "status": [wide.status, narrow.status],
                     "mc_p05": c["eta_mc"]["p05"], "mc_p95": c["eta_mc"]["p95"]})
    return rows


def grillage_rows(gr: Grillage, k: int, etas) -> list:
    rows = []
    for eta in etas:
        ind = gr.indicator(float(eta), k)
        rows.append({"eta": float(eta), "w_mm": ind["w_pred_mm"],
                     "deflection_truck_mm": ind["deflection_truck_mm"]["damaged"][k],
                     "deflection_truck_intact_mm": ind["deflection_truck_mm"]["intact"][k],
                     "EI_ratio_midspan": midspan_ratio(ind["EI_ratio_Gk"]), "EI_ratio_min": min(ind["EI_ratio_Gk"]),
                     "converged": ind["converged"]})
    return rows


def damaged_length_study(eta: float) -> dict:
    by_mesh = {"tied": [], "fixed_L5": [], "distributed": [], "intact": []}
    for n in MESH_NS:
        g = Girder1D(n_el=n)
        L = g.inp.span
        by_mesh["tied"].append(g.indicator(eta, damaged=g.damaged_zone(2 * L / n))["deflection_truck_mm"])
        by_mesh["fixed_L5"].append(g.indicator(eta, damaged=g.damaged_zone(L / 5))["deflection_truck_mm"])
        by_mesh["distributed"].append(g.indicator(eta)["deflection_truck_mm"])
        by_mesh["intact"].append(g.indicator(1.0)["deflection_truck_mm"])
    fine = g                                           # the last mesh, N = 160
    Le = fine.inp.span / fine.n_el
    lengths = [2 * Le * m for m in LENGTH_STEPS]
    by_length = [fine.indicator(eta, damaged=fine.damaged_zone(ell))["deflection_truck_mm"] for ell in lengths]
    dist = fine.indicator(eta)
    if abs(by_length[-1] - dist["deflection_truck_mm"]) > 1e-9 * dist["deflection_truck_mm"]:
        raise RuntimeError("damage over the whole span must reproduce the distributed model")
    cracked = np.array(dist["EI_ratio"]) < 1.0 - 1e-9
    x_mid = np.array(dist["x_mid"])
    return {"eta": eta, "N": list(MESH_NS),
            "tied_length_m": [2 * fine.inp.span / n for n in MESH_NS],
            "fixed_length_m": fine.inp.span / 5,
            "deflection_truck_mm": by_mesh,
            "N_fine": fine.n_el, "lengths_m": lengths, "deflection_truck_by_length_mm": by_length,
            "cracked_zone_service_m": float(cracked.sum() * Le),
            "cracked_zone_x_m": [float(x_mid[cracked].min() - Le / 2), float(x_mid[cracked].max() + Le / 2)],
            "visible_crack_zone_at_rest_m": dist["visible_crack_zone_m"],
            "deflection_limit_mm": dist["deflection_limit_mm"]}


def hand_calc() -> dict:
    """Stage-1 girder, N = 10, intact composite EI, stage-B dead load only (wearing surface and
    barriers): one element's [k] and {f}, node 2 assembled, partition, reactions, midspan."""
    g = Girder1D(n_el=10)
    L, n = g.inp.span, g.n_el
    Le, EI, q = L / n, g.EI_I, g.q_b                   # N·m², N/m acting downward
    k = beam_k(EI, Le)
    f = udl_equivalent(-q, Le)
    K = np.zeros((2 * (n + 1), 2 * (n + 1)))
    for e in range(n):
        idx = np.arange(2 * e, 2 * e + 4)
        K[np.ix_(idx, idx)] += k
    sol = solve(g.model_dead, EI)
    mid = n // 2
    v_fe, v_exact = -sol.D[2 * mid], 5 * q * L**4 / (384 * EI)
    if abs(v_fe - v_exact) > 1e-9 * v_exact:
        raise RuntimeError("nodal deflection must be exact for Hermite elements")
    R = [float(sol.R[0]), float(sol.R[2 * n])]
    if abs(sum(R) - q * L) > 1e-9 * q * L:
        raise RuntimeError("reactions must balance the load")
    kN = 1e3
    return {"span_m": L, "n_el": n, "Le_m": Le, "EI_kNm2": EI / kN, "q_kN_m": q / kN,
            "k_e_kN_m": (k / kN).tolist(),
            "f_e_kN": (f / kN).tolist(),
            "element_dofs_1based": {"element 1": [1, 2, 3, 4], "element 5": [9, 10, 11, 12],
                                    "element 6": [11, 12, 13, 14]},
            "node2_block": (K[2:4, 2:4] / kN).tolist(),
            "restrained_dofs_1based": [1, 2 * n + 1], "free_dofs": 2 * (n + 1) - 2,
            "reactions_kN": [r / kN for r in R],
            "theta_1_rad": float(sol.D[1]), "theta_1_exact_rad": float(-q * L**3 / (24 * EI)),
            "v_mid_fe_mm": v_fe * 1e3, "v_mid_exact_mm": v_exact * 1e3}


def main() -> None:
    t0 = time.perf_counter()
    res: dict = {"hand_calc": hand_calc()}
    out = ROOT / "results"
    s1 = json.loads((out / "stage1.json").read_text(encoding="utf-8"))
    s2 = json.loads((out / "stage2.json").read_text(encoding="utf-8"))

    g = Girder1D()
    eta_dec1 = eta_decompression(g.w_pred)
    base = g.indicator(1.0)["deflection_truck_mm"]
    res["stage1"] = {"n_el": g.n_el, "eta_decompression": eta_dec1,
                     "deflection_limit_mm": g.indicator(1.0)["deflection_limit_mm"],
                     "rows": stage1_rows(g, [1.0, *np.linspace(eta_dec1, ETA_MIN, N_ETA_STAGE1)]),
                     "bands": [stage1_case_band(g, c, base) for c in s1["cases"]
                               if c["status"] == "measured" and c["gsd_mm"] == BAND_GSD_MM]}
    ends = end_point_rows(g, [c for c in s1["cases"] if c["status"] == "measured"])
    for r in ends:
        # same girder as run_all.py; with sigma = u/2 the 5-95 % band (±1.645 sigma) lies inside ±u
        if abs(estimate_eta(g.w_pred, r["w_mm"] * 1e-3).eta - r["eta_point"]) > 1e-6:
            raise RuntimeError(f"{r['name']}: eta at w differs from stage1.json")
        if r["status"] != ["ok", "ok"] or not r["eta_at_w_plus_u"] <= r["mc_p05"] <= r["mc_p95"] <= r["eta_at_w_minus_u"]:
            raise RuntimeError(f"{r['name']}: the end points must bracket the Monte Carlo band")
    res["stage1"]["end_points"] = ends
    gr = Grillage(n_el=N_EL_GRID)
    res["grillage"] = {}
    for name, k in GIRDERS.items():
        eta_dec = eta_decompression(gr.w_pred_fn(k))
        ref = next(c for c in s2[name[:2]]["cases"] if c["name"] == s2[name[:2]]["reference_case"])
        res["grillage"][name] = {"k": k, "eta_decompression": eta_dec,
                                 "rows": grillage_rows(gr, k, [1.0, *np.linspace(eta_dec, ETA_MIN, N_ETA_GRID)]),
                                 "band_from_stage2": {"name": ref["name"], "w": ref["w_meas_mm"], "u": ref["u95_mm"],
                                                      "deflection_truck_mm": {
                                                          "lo": s2[name[:2]]["deflection_truck_band_mm"]["p05"],
                                                          "mid": s2[name[:2]]["deflection_truck_band_mm"]["p50"],
                                                          "hi": s2[name[:2]]["deflection_truck_band_mm"]["p95"]}}}
    res["damaged_length"] = damaged_length_study(ETA_REF)
    res["timing_s"] = time.perf_counter() - t0

    (out / "proposal_checks.json").write_text(json.dumps(res, indent=2, ensure_ascii=False, allow_nan=False),
                                              encoding="utf-8")
    make_figures(res, out / "figures")
    h, d = res["hand_calc"], res["damaged_length"]
    print("hand calc: v_mid", round(h["v_mid_fe_mm"], 4), "mm, reactions", [round(r, 2) for r in h["reactions_kN"]])
    print("damaged length (truck mm):", {k: [round(v, 2) for v in vs] for k, vs in d["deflection_truck_mm"].items()})
    print("cracked zone in service:", d["cracked_zone_service_m"], "m", d["cracked_zone_x_m"])
    print("timing", round(res["timing_s"], 1), "s")


def make_figures(res: dict, figs: Path) -> None:
    lim = res["stage1"]["deflection_limit_mm"]
    models = {"Stage-1 single girder": res["stage1"]["rows"],
              **{f"grillage {name}": v["rows"] for name, v in res["grillage"].items()}}
    gates = {f"{GATE_PX:g} px at GSD {s} mm/px": GATE_PX * s for s in GSDS_MM}

    def band(curve, b, key):
        return {"curve": curve, "w": b["w"], "u": b["u"], **b[key]}

    s1_label = "Stage-1 single girder"
    bands = [band(s1_label, b, "deflection_truck_mm") for b in res["stage1"]["bands"]]
    bands += [band(f"grillage {name}", v["band_from_stage2"], "deflection_truck_mm") for name, v in res["grillage"].items()]
    F.fig_width_response({m: ([r["w_mm"] for r in rows], [r["deflection_truck_mm"] for r in rows])
                          for m, rows in models.items()},
                         figs / "fig15_width_deflection.png", ylabel="truck-only midspan deflection (mm)",
                         hlines={f"L/800 = {lim:.1f} mm": lim}, vlines=gates, bands=bands, vline_label_y=0.36)
    s1 = res["stage1"]["rows"]
    curves = {"Stage-1 girder: midspan element EI / EI_I": ([r["w_mm"] for r in s1], [r["EI_ratio_midspan"] for r in s1]),
              "Stage-1 girder: equivalent uniform α = δ(η = 1) / δ(η)": ([r["w_mm"] for r in s1],
                                                                       [r["alpha_equivalent"] for r in s1])}
    for name, v in res["grillage"].items():
        curves[f"grillage {name}: midspan element EI / EI_I"] = ([r["w_mm"] for r in v["rows"]],
                                                                [r["EI_ratio_midspan"] for r in v["rows"]])
    labels = list(curves)
    bands = [band(labels[0], b, "EI_ratio_midspan") for b in res["stage1"]["bands"]]
    bands += [band(labels[1], b, "alpha_equivalent") for b in res["stage1"]["bands"]]
    F.fig_width_response(curves, figs / "fig16_width_stiffness.png", ylabel="stiffness ratio (service, dead + truck)",
                         vlines=gates, bands=bands)
    d = res["damaged_length"]
    by = d["deflection_truck_mm"]
    F.fig_damaged_length(d["N"], {"one zone of 2 elements at midspan (ℓ = 2L/N)": by["tied"],
                                  f"fixed ℓ = L/5 = {d['fixed_length_m']:.2f} m": by["fixed_L5"],
                                  "distributed: every element from its own moment": by["distributed"]},
                         d["lengths_m"], d["deflection_truck_by_length_mm"], figs / "fig17_damaged_length.png",
                         hlines={f"intact {by['intact'][-1]:.1f} mm": by["intact"][-1],
                                 f"L/800 = {d['deflection_limit_mm']:.1f} mm": d["deflection_limit_mm"]},
                         marks={f"cracked zone {d['cracked_zone_service_m']:.1f} m": d["cracked_zone_service_m"],
                                "proposal, N = 10: 2.45 m": d["tied_length_m"][1]})


if __name__ == "__main__":
    if "--figures-only" in sys.argv:                   # redraw from the written JSON, no analysis
        make_figures(json.loads((ROOT / "results" / "proposal_checks.json").read_text(encoding="utf-8")),
                     ROOT / "results" / "figures")
    else:
        main()
