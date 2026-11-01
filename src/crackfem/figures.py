"""Report figures. Every function takes plain arrays/dicts (from results JSON) and a path."""
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402


def _save(fig, path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _right_label(ax, y, text) -> None:
    """Small grey label just above a horizontal line, at the right edge of the axes."""
    ax.annotate(text, (1.0, y), xycoords=("axes fraction", "data"), xytext=(-3, 2),
                textcoords="offset points", ha="right", va="bottom", fontsize=7, color="0.3")


def fig_w_eta(eta, w_mm, cases, path, eta_dec=None, gate_widths_mm=None) -> None:
    """Forward model w(eta), with each measured case drawn at its measured width: eta_true (x) and
    the Monte Carlo 5-95 % band of the eta estimate. Not a validation plot: the point estimate lies
    on the curve by construction, so it is not drawn."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(eta, w_mm, color="k", label="predicted crack width at rest, w(η)")
    if eta_dec is not None:
        ax.axvline(eta_dec, color="0.4", ls="--", lw=1, label=f"η_dec = {eta_dec:.3f} (no crack predicted above)")
    for text, w in (gate_widths_mm or {}).items():
        ax.axhline(w, color="0.5", ls=":", lw=1)
        _right_label(ax, w, text)
    shown = [(i, c) for i, c in enumerate(cases) if c["status"] == "measured" and c.get("eta_mc")]
    for c in {(c["w_meas_mm"], c["u95_mm"]): c for _, c in shown}.values():   # one band per distinct (w, u)
        ax.axhspan(c["w_meas_mm"] - c["u95_mm"], c["w_meas_mm"] + c["u95_mm"], color="0.5", alpha=0.08, lw=0)
    # widest Monte Carlo band first, so the narrower ones drawn on top stay visible at a shared width
    for rank, (i, c) in enumerate(sorted(shown, key=lambda t: t[1]["eta_mc"]["p05"] - t[1]["eta_mc"]["p95"])):
        m, col = c["eta_mc"], f"C{i % 10}"
        ax.errorbar([m["p50"]], [c["w_meas_mm"]], xerr=[[m["p50"] - m["p05"]], [m["p95"] - m["p50"]]],
                    fmt="none", ecolor=col, elinewidth=2.5 - 0.3 * min(rank, 4), capsize=4,
                    label=f"{c['name']}: MC 5–95 % of η (GSD {c['gsd_mm']} mm/px)")
        ax.plot([c["eta_true"]], [c["w_meas_mm"]], "x", color=col, mew=2, ms=max(10 - 1.5 * rank, 4), zorder=5)
    handles, labels = ax.get_legend_handles_labels()
    if shown:
        handles += [Patch(color="0.5", alpha=0.25), Line2D([], [], color="k", marker="x", ls="none", mew=2)]
        labels += ["measured width ± u (u = 2 GSD)", "η true, at the measured width"]
    ax.legend(handles, labels, fontsize=6.5, loc="upper right")
    ax.set_xlabel("remaining prestress fraction η")
    ax.set_ylabel("crack width w (mm)")
    _save(fig, path)


def fig_mkappa(M_kNm, curves, path, title=None, break_at=None, hlines=None) -> None:
    """break_at: {curve label: moment (kN·m)} where that curve jumps (EC2 zeta at Mcr); the line is
    broken there instead of drawing the jump as a slope. hlines: {label: moment (kN·m)}."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for label, k in curves.items():
        x, y = [1e3 * v for v in k], [float(m) for m in M_kNm]
        mb = (break_at or {}).get(label)
        if mb is not None:
            for j in range(len(y) - 1):
                if y[j] <= mb < y[j + 1]:
                    x.insert(j + 1, math.nan)
                    y.insert(j + 1, math.nan)
                    break
        ax.plot(x, y, label=label)
    for (label, m), ls in zip((hlines or {}).items(), (":", "-.", "--")):
        ax.axhline(m, color="0.4", ls=ls, lw=1, label=label)
    ax.set_xlabel("curvature from the end of stage A (1/km)")
    ax.set_ylabel("total moment (kN·m)")
    if title:
        ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="lower right")
    _save(fig, path)


def fig_ei_ratio(x, ratios, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for label, r in ratios.items():
        ax.step(x, r, where="mid", label=label)
    ax.set_xlabel("x along the span (m)")
    ax.set_ylabel("secant EI / uncracked EI\n(EC2, stage-B dead + truck)")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_deflection(x, shapes, path, title=None) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    for label, v in shapes.items():
        ax.plot(x, v, label=label)
    ax.invert_yaxis()
    ax.set_xlabel("x along the span (m)")
    ax.set_ylabel("deflection (mm, down), EC2\nstage-B dead + truck")
    if title:
        ax.set_title(title, fontsize=8.5)
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_eta_hist(samples, eta_true, path, title=None, p05=None, p95=None) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.hist(samples, bins=40, color="0.6")
    ax.axvline(eta_true, color="k", linestyle="--", label="η true")
    if p05 is not None and p95 is not None:
        ax.axvline(p05, color="C0", linestyle=":", lw=1.5, label="MC 5 % and 95 %")
        ax.axvline(p95, color="C0", linestyle=":", lw=1.5)
    ax.set_xlabel("η estimate (Monte Carlo on w)")
    ax.set_ylabel("count")
    if title:
        ax.set_title(title, fontsize=8.5)
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_mesh(Ns, series, path, note=None, ylabel=None) -> None:
    """Relative change of each series against its finest mesh, |d_N - d_ref| / |d_ref|, on log axes.
    Points that equal the reference exactly cannot be drawn on a log axis and are left out.
    ylabel: replaces the default (Stage-1 midspan-deflection) axis label."""
    fig, ax = plt.subplots(figsize=(6, 3.5))
    n_ref = Ns[-1]
    for label, v in series.items():
        ref = v[-1]
        pts = [(n, abs(d - ref) / abs(ref)) for n, d in zip(Ns[:-1], v[:-1])]
        pts = [(n, r) for n, r in pts if r > 0.0]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(list(Ns[:-1]))
    ax.set_xticklabels([str(n) for n in Ns[:-1]])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("number of elements N")
    ax.set_ylabel(ylabel if ylabel else f"|δ_N − δ_{n_ref}| / δ_{n_ref}\nmidspan, stage-B dead + truck")
    if note:
        ax.set_title(note, fontsize=8)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7)
    _save(fig, path)


def fig_cases(cases, path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    vals = [c["eta_true"] for c in cases]
    vals += [c["eta_mc"][k] for c in cases if c["status"] == "measured" and c["eta_mc"] for k in ("p05", "p95")]
    vals += [c["eta_lower_bound_from_gate"] for c in cases if c.get("eta_lower_bound_from_gate") is not None]
    arrow = 0.08 * ((max(vals) - min(vals)) or 0.5)
    h = {}
    for i, c in enumerate(cases):
        h["η true"] = ax.plot([i], [c["eta_true"]], "kx", ms=7, mew=1.5, zorder=5)[0]
        if c["status"] == "measured" and c["eta_mc"]:
            m = c["eta_mc"]
            h["MC median (5–95 %)"] = ax.errorbar([i], [m["p50"]], yerr=[[m["p50"] - m["p05"]], [m["p95"] - m["p50"]]],
                                                  fmt="o", color="C0", capsize=3)
            continue
        lb = c.get("eta_lower_bound_from_gate")
        if lb is None:
            # "unresolved": an unreadable crack may be wider than the gate, so no bound exists either
            label = "unresolved,\nno bound" if c["status"] == "unresolved" else "refused,\nno bound"
            ax.annotate(label, (i, c["eta_true"]), textcoords="offset points", xytext=(0, 8),
                        ha="center", fontsize=6.5)
        else:
            h["refused: η ≥ bound at the 3 px gate"] = ax.plot([i], [lb], "_", color="C3", ms=12, mew=2)[0]
            ax.annotate("", xy=(i, lb + arrow), xytext=(i, lb),
                        arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.2, shrinkA=0, shrinkB=0))
    ax.margins(y=0.15)                     # room for the labels and arrows inside the frame
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([f"{c['name']}" for c in cases], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("η")
    ax.legend(list(h.values()), list(h.keys()), fontsize=7, loc="upper right")
    _save(fig, path)


def fig_tornado(rows, path, labels=None, title=None) -> None:
    """Each row spans min(low, base, high) to max: the upper half-bar runs from the base value to the
    metric at the low end of the tested range, the lower half-bar to the metric at the high end
    (both halves can lie on one side of the base when the effect is not monotone)."""
    fig, ax = plt.subplots(figsize=(6.5, 1.0 + 0.5 * len(rows)))
    base = rows[0]["base"]
    for i, r in enumerate(reversed(rows)):
        ax.barh(i + 0.19, r["metric_low"] - base, left=base, height=0.38, color="0.35")
        ax.barh(i - 0.19, r["metric_high"] - base, left=base, height=0.38, color="0.78")
    ax.axvline(base, color="k")
    ax.use_sticky_edges = False            # bars' ends must not sit on the frame
    ax.margins(x=0.05)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([(labels or {}).get(r["name"], r["name"]) for r in reversed(rows)], fontsize=8)
    ax.set_xlabel("η estimate")
    ax.legend([Patch(color="0.35"), Patch(color="0.78"), Line2D([], [], color="k")],
              ["low end of the tested range", "high end of the tested range", f"base η = {base:.3f}"],
              fontsize=7, loc="lower right")
    if title:
        ax.set_title(title, fontsize=9)
    _save(fig, path)


# --- Stage 2 (grillage) ------------------------------------------------------------------
ETA_GEN_LABEL = "η used to generate the crack (1-D model)"


def fig_w_eta_models(eta, curves_mm, path, eta_dec=None) -> None:
    """Predicted crack width at rest w(eta) of several models (e.g. Stage-1 girder, grillage G1,
    grillage G3); eta_dec: {label: plateau edge} drawn as dashed vertical lines."""
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    for i, (label, w) in enumerate(curves_mm.items()):
        ax.plot(eta, w, color=f"C{i}", label=label)
        if eta_dec and label in eta_dec:
            ax.axvline(eta_dec[label], color=f"C{i}", ls="--", lw=1)
    ax.set_xlabel("remaining prestress fraction η of the cracked girder")
    ax.set_ylabel("crack width at rest w (mm)")
    handles, names = ax.get_legend_handles_labels()
    if eta_dec:
        handles.append(Line2D([], [], color="0.4", ls="--", lw=1))
        names.append("dashed: η_dec, no crack at rest above it")
    ax.legend(handles, names, fontsize=7)
    _save(fig, path)


def fig_eta_bands(rows, path, title=None) -> None:
    """rows: [{"name", "eta_gen_1d", "bands": {model label: {"p05", "p50", "p95"} or None}}]: the
    Monte Carlo 5-50-95 % of eta per measured case and model, side by side. eta_gen_1d is the eta
    the synthetic crack was generated with (Stage-1 single-girder model), not a truth for the
    grillage."""
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    labels = sorted({m for r in rows for m in r["bands"]}, key=lambda s: (s != "Stage 1", s))
    width = 0.8 / max(len(labels), 1)
    handles = {}
    for i, r in enumerate(rows):
        handles[ETA_GEN_LABEL] = ax.plot([i - 0.45, i + 0.45], [r["eta_gen_1d"]] * 2, "k--", lw=1, zorder=1)[0]
        for j, m in enumerate(labels):
            b = r["bands"].get(m)
            if not b:
                continue
            xj = i - 0.4 + width * (j + 0.5)
            handles[m] = ax.errorbar([xj], [b["p50"]], yerr=[[b["p50"] - b["p05"]], [b["p95"] - b["p50"]]],
                                     fmt="o", color=f"C{j}", capsize=3, ms=4)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([r["name"] for r in rows], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("η (Monte Carlo 5–50–95 %)")
    if title:
        ax.set_title(title, fontsize=9)
    ax.legend(list(handles.values()), list(handles.keys()), fontsize=7, loc="upper right")
    _save(fig, path)


def fig_girder_bars(girders, series, path, ylabel, title=None, hline=None) -> None:
    """Grouped bars per girder (e.g. intact vs damaged midspan moment or deflection); hline:
    optional (value, label) drawn as a dashed horizontal line (e.g. a deflection limit)."""
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    n = len(series)
    width = 0.8 / max(n, 1)
    for j, (label, vals) in enumerate(series.items()):
        xs = [i - 0.4 + width * (j + 0.5) for i in range(len(girders))]
        ax.bar(xs, vals, width=width, color=f"C{j}", label=label)
    ax.axhline(0.0, color="k", lw=0.8)
    if hline:
        ax.axhline(hline[0], color="0.3", ls="--", lw=1, label=hline[1])
        ax.set_ylim(top=1.15 * max(hline[0], max(max(v) for v in series.values())))
    ax.set_xticks(range(len(girders)))
    ax.set_xticklabels(girders)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_history(series, path, title=None) -> None:
    """Relative moment change after each secant update (log scale), one line per analysis."""
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    for label, h in series.items():
        pts = [(i + 1, v) for i, v in enumerate(h) if v > 0.0]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=3, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("stiffness update")
    ax.set_ylabel("max |ΔM| / max |M|")
    if title:
        ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
    _save(fig, path)


NO_CRACK_MM = 1e-3


def fig_width_response(curves, path, ylabel, hlines=None, vlines=None, title=None, bands=None,
                       vline_label_y: float = 1.0) -> None:
    """A response against the measured crack width at rest, one curve per model.
    curves: {label: (w_mm, values)}; leading points with w below NO_CRACK_MM (no crack at rest:
    the blind zone; at eta_dec the width is a few nanometres) are joined by a dotted line. hlines: {label: value}; vlines: {label: width (mm)}, e.g. the smallest width
    the pipeline measures at a GSD (3 px), labelled at vline_label_y (axes fraction; 1 = top).
    bands: measured cases carried through, [{"curve": label, "w": w, "u": u, "lo": p05, "mid": p50,
    "hi": p95}], drawn as w ± u across and the Monte Carlo 5-95 % band of the response up and
    down, in the colour of that curve."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    colour = {label: f"C{i}" for i, label in enumerate(curves)}
    def n_blind(w):
        return next((i for i, x in enumerate(w) if x >= NO_CRACK_MM), len(w))

    for label, (w, v) in curves.items():
        n0 = n_blind(w)
        if n0 > 1:                                   # blind zone: several eta give w = 0
            ax.plot(w[:n0], v[:n0], "o:", ms=3, color=colour[label])
        ax.plot(w[max(n0 - 1, 0):], v[max(n0 - 1, 0):], "o-", ms=3, color=colour[label], label=label)
    if any(n_blind(w) > 1 for w, _ in curves.values()):
        ax.plot([], [], "o:", ms=3, color="0.3", label="w = 0: no crack at rest, from η = 1 to η_dec (blind zone)")
    for b in bands or []:
        ax.errorbar([b["w"]], [b["mid"]], xerr=[[b["u"]], [b["u"]]],
                    yerr=[[b["mid"] - b["lo"]], [b["hi"] - b["mid"]]],
                    fmt="s", ms=4, color=colour[b["curve"]], ecolor=colour[b["curve"]], capsize=3, zorder=5)
    if bands:
        ax.plot([], [], "s", color="0.3", ms=4, label="measured case: w ± u, MC 5–95 % of the response")
    for text, y in (hlines or {}).items():
        ax.axhline(y, color="0.3", ls="--", lw=1)
        _right_label(ax, y, text)
    top = vline_label_y >= 1.0
    for text, x in (vlines or {}).items():
        ax.axvline(x, color="0.5", ls=":", lw=1)
        ax.annotate(text, (x, vline_label_y), xycoords=("data", "axes fraction"), xytext=(2, -3 if top else 0),
                    textcoords="offset points", rotation=90, ha="left", va="top" if top else "bottom",
                    fontsize=6.5, color="0.3")
    ax.set_xlim(left=-0.01)
    ax.set_xlabel("measured crack width at rest w (mm)")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=8)
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_damaged_length(Ns, by_mesh, lengths_m, by_length, path, hlines=None, marks=None) -> None:
    """Left: a response against the number of elements N for damage zones defined in different
    ways ({label: values per N}). Right: the same response against a fixed damaged length on a
    fine mesh. hlines: {label: value} on both panels; marks: {label: length (m)} on the right."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.0, 3.8), sharey=True)
    for i, (label, v) in enumerate(by_mesh.items()):
        a.plot(Ns, v, "o-", color=f"C{i}", label=label)
    a.set_xscale("log")
    a.set_xticks(list(Ns))
    a.set_xticklabels([str(n) for n in Ns])
    a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel("number of elements N")
    b.plot(lengths_m, by_length, "o-", color="k", ms=3, label="damage confined to length ℓ around midspan")
    for text, x in (marks or {}).items():
        b.axvline(x, color="0.5", ls=":", lw=1)
        b.annotate(text, (x, 0.45), xycoords=("data", "axes fraction"), xytext=(2, 0),
                   textcoords="offset points", rotation=90, ha="left", va="bottom", fontsize=6.5, color="0.3")
    b.set_xlabel("damaged length ℓ (m)")
    for ax in (a, b):
        for text, y in (hlines or {}).items():
            ax.axhline(y, color="0.3", ls="--", lw=1)
            _right_label(ax, y, text)
    a.legend(fontsize=6.5, loc="center left", bbox_to_anchor=(0.0, 0.68))
    b.legend(fontsize=6.5, loc="lower right", bbox_to_anchor=(1.0, 0.12))
    a.set_ylabel("truck-only midspan deflection (mm)")
    _save(fig, path)
