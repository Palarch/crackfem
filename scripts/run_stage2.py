"""Regenerate every Stage-2 (grillage) result and figure from committed inputs.

Reads data/claimbridge/*_spike.json (thesis-pipeline output) and results/stage1.json (Stage-1
bands, read only) and writes results/stage2.json + results/figures/fig9...fig14. Stage-1 outputs
are not touched. Does not need the thesis pipeline; needs PyNiteFEA (requirements.txt).

The synthetic crack widths were generated with the Stage-1 single-girder model, so the eta they
were made with is reported as `eta_gen_1d` ("eta used to generate the crack (1-D model)"): it is
a Stage-1 value, not a truth for the grillage."""
import json
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
from crackfem.gridsolve import moment_at, solve  # noqa: E402
from crackfem.grillage import Grillage  # noqa: E402
from crackfem.inverse import estimate_eta  # noqa: E402
from crackfem.params import Inputs  # noqa: E402
from crackfem.uq import monte_carlo_eta  # noqa: E402
from crackfem.wtable import build_table, max_rel_error  # noqa: E402

N_EL = 20
MC_N = 2000
GRID = np.linspace(0.02, 1.0, 99)
MESH_NS = (10, 20, 40, 80)
GJ_SENSITIVITY = 0.5
CRACKED = {"G1": 0, "G3": 2}
GATE_PX = 3.0
QUANTILES = ("p05", "p50", "p95")
# G1 model sensitivity at its reference eta: (diaphragm, lane_rule); ("web", "lanes") is the default
MODEL_VARIANTS = (("web", "lanes"), ("none", "lanes"), ("composite", "lanes"), ("web", "clearance"))


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return jsonable(o.item())
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def run_case(table, w_exact, rec: dict) -> dict:
    """One synthetic case on one grillage model: exact point estimate, Monte Carlo on the table,
    the table's percentiles confirmed by the exact inverse; refusals give the nominal gate bound.
    eta_gen_1d: the eta the synthetic crack was generated with (Stage-1 single-girder model)."""
    c = rec["case"]
    out = {"name": c["name"], "gsd_mm": c["gsd_mm"], "eta_gen_1d": c["eta_true"]}
    gov = governing_width(rec["parsed"])
    if gov is None:
        if not refused_below_resolution(rec["parsed"]):
            return {**out, "status": "unresolved", "eta_mc": None, "eta_lower_bound_from_gate": None}
        gate = estimate_eta(w_exact, GATE_PX * c["gsd_mm"] * 1e-3)
        return {**out, "status": "refused", "eta_mc": None, "gate_inverse_status": gate.status,
                "eta_lower_bound_from_gate": gate.eta if gate.status == "ok" else None}
    w, u = gov["max_width_mm"] * 1e-3, gov["uncertainty_mm"] * 1e-3
    est = estimate_eta(w_exact, w)
    mc = monte_carlo_eta(table, w, u, n=MC_N)
    s = mc.summary()
    exact = {}
    if s:
        for q in QUANTILES:
            e = estimate_eta(w_exact, table(s[q]))
            exact[q] = e.eta if e.status == "ok" else None
    return {**out, "status": "measured", "w_meas_mm": gov["max_width_mm"], "u95_mm": gov["uncertainty_mm"],
            "eta_point": est.eta, "eta_status": est.status, "eta_point_evaluations": est.evaluations,
            "eta_mc": s, "eta_mc_exact_check": exact,
            "eta_mc_exact_max_abs_diff": (max(abs(exact[q] - s[q]) for q in QUANTILES)
                                          if s and all(exact[q] is not None for q in QUANTILES) else None)}


def pynite_crosscheck(g: Grillage) -> dict:
    from crackfem.pynite_check import pynite_deflections, pynite_grid, pynite_moment

    model = g.service_model(0)
    ours = solve(model, g.EI_I, g.GJ0)
    fm = pynite_grid(model, g.EI_I, g.GJ0)
    w_ours, w_pyn = ours.D[0::3], pynite_deflections(fm, model)
    d = g.deck
    right_of_mid = [int(d.girder_members[k, d.mid]) for k in range(d.n_girders)]   # members starting at midspan
    m_ours = np.array([moment_at(model, ours, m, 0.0) for m in right_of_mid])
    m_pyn = np.array([pynite_moment(fm, m, 0.0) for m in right_of_mid])
    sign = float(np.sign(m_pyn[2] / m_ours[2]))                  # PyNite's own Mz sign convention
    return {"model": "linear intact grillage, stage-B dead + truck at G1's critical placement",
            "max_abs_w_diff_mm": float(np.max(np.abs(w_ours - w_pyn)) * 1e3),
            "max_w_mm": float(np.max(np.abs(w_ours)) * 1e3),
            "pynite_Mz_sign": sign,
            "max_rel_M_diff": float(np.max(np.abs(sign * m_pyn - m_ours)) / np.max(np.abs(m_ours))),
            "M_mid_right_ours_kNm": m_ours / 1e3, "M_mid_right_pynite_kNm": sign * m_pyn / 1e3}


def dead_moments(ind: dict) -> dict:
    """Stage-B dead-load (at-rest) midspan moments: service minus truck part."""
    return {s: [a - b for a, b in zip(ind["M_B_service_mid_kNm"][s], ind["M_truck_mid_kNm"][s])]
            for s in ("intact", "damaged")}


def totals(g: Grillage, ind: dict) -> dict:
    """Total midspan moments (stage A from statics + stage B from the grillage) of every girder,
    at rest (stage-B dead load) and in service (+ truck), intact and damaged. Stage A is the same
    intact or damaged: it stays on each girder (statics on the bare girder)."""
    M_A = [float(gd.M_A_nd[g.deck.mid]) / 1e3 for gd in g.girders]
    dead = dead_moments(ind)
    return {"M_A_mid_kNm": M_A,
            "M_total_mid_at_rest_kNm": {s: [a + b for a, b in zip(M_A, dead[s])] for s in ("intact", "damaged")},
            "M_total_mid_service_kNm": {s: [a + b for a, b in zip(M_A, ind["M_B_service_mid_kNm"][s])]
                                        for s in ("intact", "damaged")}}


def at_eta(g: Grillage, eta: float, k: int) -> dict:
    """Indicator of girder k cracked at a given eta, for a like-for-like comparison with Stage 1."""
    ind = g.indicator(eta, k)
    return {"eta": eta, "w_pred_mm": ind["w_pred_mm"],
            "deflection_truck_mm": ind["deflection_truck_mm"],
            "deflection_service_mm": ind["deflection_service_mm"],
            "M_B_service_mid_kNm": ind["M_B_service_mid_kNm"],
            "shed_percent_of_Gk": ind["shed_percent_of_Gk"],
            "Gk_loses_percent": -ind["shed_percent_of_Gk"][k],
            "cracked_elements": ind["cracked_elements"], "precracked_count": ind["precracked_count"],
            "converged": ind["converged"], **totals(g, ind)}


def model_variant(gv: Grillage, eta: float, ref_rec: dict, strands: str) -> dict:
    """G1 cracked at a fixed eta on one model variant, and the reference case inverted on it."""
    ind = gv.indicator(eta, 0)
    f = gv.w_pred_fn(0)
    table = build_table(f, GRID)
    case = run_case(table, f, ref_rec)
    return {"diaphragm": gv.diaphragm, "lane_rule": gv.lane_rule, "strands": strands,
            "n_strands": gv.n_strands, "eta": eta,
            "M_truck_mid_critical_kNm": [gv.M_truck_mid[c, k] / 1e3 for k, c in enumerate(gv.critical)],
            "G1_loses_percent": -ind["shed_percent_of_Gk"][0],
            "G2_takes_percent_of_G1": ind["shed_percent_of_Gk"][1],
            "w_pred_mm": ind["w_pred_mm"], "cracked_elements": ind["cracked_elements"],
            "G2_cracked_elements": ind["cracked_elements"][1],
            "deflection_truck_mm": ind["deflection_truck_mm"],
            "converged": ind["converged"], "eta_decompression": table.eta_dec,
            "reference_case": case["name"], "eta_point_ref_case": case.get("eta_point"),
            "eta_mc_ref_case": case.get("eta_mc")}


def sensitivity_G1(g: Grillage, eta: float, ref_rec: dict) -> list:
    """Diaphragm none / web / composite and lane rule clearance for G1 at its reference eta, strands
    held at the default design; a variant whose own zero-tension design differs is run again
    with its own strands."""
    rows = []
    for diaphragm, lane in MODEL_VARIANTS:
        if (diaphragm, lane) == (g.diaphragm, g.lane_rule):
            gv = g
        else:
            gv = Grillage(n_el=N_EL, n_strands=g.n_strands, diaphragm=diaphragm, lane_rule=lane)
        row = model_variant(gv, eta, ref_rec, "default design")
        own = g if gv is g else Grillage(n_el=N_EL, diaphragm=diaphragm, lane_rule=lane)
        row["own_design_n_strands"] = own.n_strands
        rows.append(row)
        if own.n_strands != g.n_strands:
            rows.append(model_variant(own, eta, ref_rec, "own design"))
    return rows


def main(pynite: bool = True) -> None:
    t_all = time.perf_counter()
    res: dict = {"inputs": asdict(Inputs()), "n_el": N_EL,
                 "note": ("eta_gen_1d = the eta used to generate the synthetic crack with the Stage-1 "
                          "single-girder (1-D) model; it is not a truth for the grillage")}
    t0 = time.perf_counter()
    g = Grillage(n_el=N_EL)
    t_build = time.perf_counter() - t0
    res["model"] = {"diaphragm": g.diaphragm, "lane_rule": g.lane_rule, "n_steps": g.n_steps}
    res["design"] = {
        "n_strands": g.n_strands, "design_rounds": g.design_rounds, "slab_width_m": g.widths,
        "soffit_stress_service_MPa": [s / 1e6 for s in g.soffit_service],
        "J_girder_m4": g.J_girder, "q_stage_b_dead_kN_m": g.q_dead / 1e3, "torque_stage_b_dead_kNm_m": g.t_dead / 1e3,
        "critical_placement_wheel_lines_m": [g.placements[c] for c in g.critical],
        "n_placements": len(g.placements),
        "M_B_dead_mid_intact_kNm": g.M_B_dead_mid_I / 1e3,
        "M_truck_mid_critical_kNm": [g.M_truck_mid[c, k] / 1e3 for k, c in enumerate(g.critical)],
        "truck_scale": g.truck_scale}
    stage1 = json.loads((ROOT / "results" / "stage1.json").read_text(encoding="utf-8"))
    s1_cases = {c["name"]: c for c in stage1["cases"]}
    g1d = Girder1D()
    res["stage1_w_curve"] = {"eta": GRID, "w_mm": [g1d.w_pred(float(e)) * 1e3 for e in GRID]}
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "data" / "claimbridge").glob("*_spike.json"))]
    recs_by_name = {r["case"]["name"]: r for r in recs}

    timing = {"build_grillage_s": t_build}
    for name, k in CRACKED.items():
        f = g.w_pred_fn(k)
        t0 = time.perf_counter()
        table = build_table(f, GRID)
        timing[f"table_{name}_s"] = time.perf_counter() - t0
        cracked_knots = table.eta[table.eta < table.eta_dec]
        mids = 0.5 * (cracked_knots[:-1] + cracked_knots[1:])
        t0 = time.perf_counter()
        cases = [run_case(table, f, r) for r in recs]
        timing[f"cases_{name}_s"] = time.perf_counter() - t0
        measured = [c for c in cases if c["status"] == "measured" and c["eta_mc"]]
        # reference case as in Stage 1: the finest GSD, and among those the least damage (the largest
        # generating eta, i.e. the first detectable case); the key is unique for the committed cases
        ref = min(measured, key=lambda c: (c["gsd_mm"], -c["eta_gen_1d"]))
        eta_ref = ref["eta_mc"]["p50"]
        band = {q: g.indicator(ref["eta_mc"][q], k) for q in QUANTILES}
        from_eta = {"p05": "p95", "p50": "p50", "p95": "p05"}       # deflection falls as eta rises
        dmg = band["p50"]
        sens = g.indicator(eta_ref, k, gj_factor=GJ_SENSITIVITY)
        res[name] = {
            "k": k, "eta_decompression": table.eta_dec,
            "w_curve": {"eta": [float(e) for e in GRID], "w_mm": [table(float(e)) * 1e3 for e in GRID]},
            "table": {"knots": table.eta.size, "max_rel_error_at_midpoints": max_rel_error(table, f, mids)},
            "cases": cases, "reference_case": ref["name"], "eta_ref": eta_ref,
            "eta_ref_stage1": s1_cases[ref["name"]]["eta_mc"]["p50"],
            "indicator": dmg,
            "M_B_dead_mid_kNm": dead_moments(dmg),
            **totals(g, dmg),
            "deflection_service_band_mm": {q: band[v]["deflection_service_mm"]["damaged"][k] for q, v in from_eta.items()},
            "deflection_truck_band_mm": {q: band[v]["deflection_truck_mm"]["damaged"][k] for q, v in from_eta.items()},
            "gj_sensitivity": {"gj_factor": GJ_SENSITIVITY,
                               "M_B_service_mid_kNm": sens["M_B_service_mid_kNm"]["damaged"],
                               "deflection_service_mm": sens["deflection_service_mm"]["damaged"],
                               "w_pred_mm": sens["w_pred_mm"]},
        }
    eta1 = float(stage1["eta_ref"])
    s1i = stage1["indicator_ref"]
    t0 = time.perf_counter()
    res["at_stage1_eta"] = {
        "eta": eta1,
        "note": "every model at the same eta (Stage-1 eta_ref): like-for-like deflections and moments",
        "stage1": {"deflection_truck_mm": s1i["deflection_truck_mm"],
                   "deflection_service_mm": s1i["deflection_service_mm"],
                   "M_total_mid_at_rest_kNm": s1i["M_dead_mid_kNm"], "w_pred_mm": s1i["w_pred_mm"]},
        **{name: at_eta(g, eta1, k) for name, k in CRACKED.items()}}
    timing["at_stage1_eta_s"] = time.perf_counter() - t0
    res["stage1_comparison"] = [
        {"name": c["name"], "eta_gen_1d": c["eta_gen_1d"],
         "bands": {"Stage 1": s1_cases[c["name"]].get("eta_mc"),
                   **{f"grillage {n}": next(x for x in res[n]["cases"] if x["name"] == c["name"])["eta_mc"] for n in CRACKED}}}
        for c in res["G1"]["cases"] if c["status"] == "measured"]

    t0 = time.perf_counter()
    res["sensitivity_G1"] = {"eta": res["G1"]["eta_ref"], "reference_case": res["G1"]["reference_case"],
                             "default_n_strands": g.n_strands,
                             "variants": sensitivity_G1(g, res["G1"]["eta_ref"],
                                                        recs_by_name[res["G1"]["reference_case"]])}
    timing["sensitivity_G1_s"] = time.perf_counter() - t0

    if pynite:
        res["pynite"] = pynite_crosscheck(g)

    t0 = time.perf_counter()
    mesh = {"N": list(MESH_NS)}
    for name, k in CRACKED.items():
        eta = res[name]["eta_ref"]
        w, dfl = [], []
        for n in MESH_NS:
            gm = g if n == N_EL else Grillage(n_el=n, n_strands=g.n_strands)
            w.append(gm.w_pred(eta, k) * 1e3)
            dfl.append(gm.indicator(eta, k)["deflection_service_mm"]["damaged"][k])
        mesh[name] = {"eta": eta, "w_pred_mm": w, "deflection_service_mm": dfl,
                      "rel_to_N80_w_pred": [abs(v - w[-1]) / abs(w[-1]) for v in w],
                      "rel_to_N80_deflection_service": [abs(v - dfl[-1]) / abs(dfl[-1]) for v in dfl]}
    res["mesh"] = mesh
    timing["mesh_s"] = time.perf_counter() - t0

    res["convergence"] = {n: {"converged": res[n]["indicator"]["converged"],
                              "iterations_last_step": res[n]["indicator"]["iterations"],
                              "passes_per_step": res[n]["indicator"]["passes_per_step"],
                              "every_table_point_converged": True}       # build_table raises otherwise
                          for n in CRACKED}

    t0 = time.perf_counter(); solve(g.service_model(0), g.EI_I, g.GJ0); timing["linear_solve_s"] = time.perf_counter() - t0
    t0 = time.perf_counter(); g.analyse(g.service_model(0), g.etas(0, res["G1"]["eta_ref"])); timing["nonlinear_service_G1_s"] = time.perf_counter() - t0
    timing["total_s"] = time.perf_counter() - t_all
    res["timing_s"] = timing

    out = ROOT / "results"
    (out / "stage2.json").write_text(json.dumps(jsonable(res), indent=2, ensure_ascii=False, allow_nan=False),
                                     encoding="utf-8")
    make_figures(json.loads((out / "stage2.json").read_text(encoding="utf-8")), stage1)

    for name in CRACKED:
        r = res[name]
        print(name, "eta_dec", round(r["eta_decompression"], 4), "eta_ref", round(r["eta_ref"], 4),
              "table err", r["table"]["max_rel_error_at_midpoints"])
        for c in r["cases"]:
            print("  ", c["name"], c["status"], c.get("eta_mc") and {q: round(c["eta_mc"][q], 4) for q in QUANTILES},
                  c.get("eta_mc_exact_max_abs_diff"), c.get("eta_lower_bound_from_gate"))
    for v in res["sensitivity_G1"]["variants"]:
        print("sens G1", v["diaphragm"], v["lane_rule"], v["strands"], v["n_strands"],
              "loses %", round(v["G1_loses_percent"], 1), "w_pred mm", round(v["w_pred_mm"], 4),
              "G2 cracked", v["G2_cracked_elements"], "eta p50", v["eta_mc_ref_case"] and round(v["eta_mc_ref_case"]["p50"], 4))
    a = res["at_stage1_eta"]
    print("at Stage-1 eta", a["eta"], "truck deflection mm: Stage 1", round(a["stage1"]["deflection_truck_mm"], 2),
          *[f"{n} {a[n]['deflection_truck_mm']['damaged'][k]:.2f} (intact {a[n]['deflection_truck_mm']['intact'][k]:.2f})"
            for n, k in CRACKED.items()])
    print("design", res["design"]["n_strands"], "timing", res["timing_s"])
    if pynite:
        print("pynite", res["pynite"]["max_abs_w_diff_mm"], res["pynite"]["max_w_mm"], res["pynite"]["max_rel_M_diff"])


def make_figures(res: dict, stage1: dict) -> None:
    """fig9-fig14 from the stage2.json content (and the Stage-1 plateau edge, read only)."""
    figs = ROOT / "results" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    s1_label = "Stage 1: 1-D interior girder, S/5.5 (generated the synthetic cracks)"
    F.fig_w_eta_models(res["stage1_w_curve"]["eta"], {s1_label: res["stage1_w_curve"]["w_mm"],
                                                      "grillage, G1 cracked": res["G1"]["w_curve"]["w_mm"],
                                                      "grillage, G3 cracked": res["G3"]["w_curve"]["w_mm"]},
                       figs / "fig9_grillage_w_eta.png",
                       eta_dec={s1_label: stage1["eta_decompression"],
                                "grillage, G1 cracked": res["G1"]["eta_decompression"],
                                "grillage, G3 cracked": res["G3"]["eta_decompression"]})
    F.fig_eta_bands(res["stage1_comparison"], figs / "fig10_eta_bands.png",
                    title="η per measured case: Stage-1 girder vs grillage (G1, G3 cracked)")
    n_girders = len(res["design"]["n_strands"])
    labels = [f"G{i + 1}" for i in range(n_girders)]
    ratio = res["inputs"]["deflection_limit_ratio"]
    for name, fig_m, fig_d in (("G1", "fig11_moment_share_G1.png", "fig12_deflection_G1.png"),
                               ("G3", "fig11_moment_share_G3.png", "fig12_deflection_G3.png")):
        ind, eta = res[name]["indicator"], res[name]["eta_ref"]
        F.fig_girder_bars(labels, {"intact": ind["M_B_service_mid_kNm"]["intact"],
                                   f"{name} at η = {eta:.3f}": ind["M_B_service_mid_kNm"]["damaged"]},
                          figs / fig_m, "stage-B midspan moment (kN·m)\ndead + truck",
                          title=(f"{name} cracked: stage-B load shed to the other girders\n"
                                 "(stage A, from statics, stays on each girder unchanged)"))
        lim = ind["deflection_limit_mm"]
        F.fig_girder_bars(labels, {"intact": ind["deflection_truck_mm"]["intact"],
                                   f"{name} at η = {eta:.3f}": ind["deflection_truck_mm"]["damaged"]},
                          figs / fig_d, "truck-only midspan deflection\n(mm, down)",
                          title=f"{name} cracked: truck deflection against L/{ratio:g}",
                          hline=(lim, f"limit L/{ratio:g} = {lim:.1f} mm"))
    F.fig_history({f"{n} service (η = {res[n]['eta_ref']:.3f})": res[n]["indicator"]["history_service"] for n in CRACKED}
                  | {f"{n} dead (at rest)": res[n]["indicator"]["history_dead"] for n in CRACKED},
                  figs / "fig13_secant_history.png",
                  title="secant iteration on the grillage (last load step, damaged analyses)")
    N = res["mesh"]["N"]
    F.fig_mesh(N, {f"w_pred at rest, {n}": res["mesh"][n]["w_pred_mm"] for n in CRACKED}
               | {f"service deflection of {n} (damaged)": res["mesh"][n]["deflection_service_mm"] for n in CRACKED},
               figs / "fig14_grillage_mesh.png", note=f"grillage mesh: relative change against N = {N[-1]}",
               ylabel=f"|v_N − v_{N[-1]}| / v_{N[-1]}\n(v = w_pred or midspan deflection)")


if __name__ == "__main__":
    if "--figures-only" in sys.argv[1:]:          # redraw fig9-fig14 from results/stage2.json
        make_figures(json.loads((ROOT / "results" / "stage2.json").read_text(encoding="utf-8")),
                     json.loads((ROOT / "results" / "stage1.json").read_text(encoding="utf-8")))
    else:
        main()
