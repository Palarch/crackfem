# Plan 2 — Stage 2: Grillage of CPK72 S01 (B4–B5)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put the Stage-1 girders into a grillage of the whole S01 deck (five girders, slab strips, diaphragms) so that the secant-stiffness loop and the η inverse run on an indeterminate FE model: a crack width on the exterior girder G1 (and, as a second case, the interior girder G3) gives η with a Monte Carlo band, and the damaged girder's stiffness loss shows up as moment shed to its neighbours and as deflections of every girder.

**Architecture:** A grid element (Logan §5.4) and a grid assembly/partition solve sit next to the beam ones and know nothing about concrete. `deck.py` turns the S01 inputs into a grillage (girder lines, slab strips, diaphragms, Hambly properties, lever-rule load mapping, truck placements). `grillage.py` builds one Stage-1 `Girder1D` per girder line (its own slab width and strand count) and reuses its staged sections, Mcr(η) and crack width unchanged; the generic `nonlinear.secant_iterate` updates every girder element's EI from its own moment and returns non-convergence as a result. `wtable.py` tabulates the expensive forward map w(η) for the Monte Carlo; `uq.monte_carlo_eta` and `inverse.estimate_eta` are reused unchanged. `pynite_check.py` solves the same linear grid in PyNite. `scripts/run_stage2.py` writes `results/stage2.json` and fig9–fig14; Stage-1 code and outputs are not changed.

**Tech Stack:** Python 3.13, numpy 2.3.5, scipy 1.17.1, sympy 1.14.0, pillow 12.0.0, pytest 8.3.5, matplotlib 3.11.2 (all pinned already, unchanged), PyNiteFEA 2.0.2 (approved dependency, installed and pinned in Task 8).

**Spec:** `docs/superpowers/specs/2026-10-07-crack-informed-grillage-design.md` (with its "As-built amendments"). **Plan 1:** `docs/superpowers/plans/2026-10-07-plan1-stage1-girder.md` (its Conventions hold here). **Scope:** B4 (grid element, deck model, nonlinear grillage, PyNite cross-check) and B5 (grillage inverse, Monte Carlo, results, figures). The report (B6) is not in this plan.

---

## Starting point and execution notes (2026-10-15)

- Work on the local branch `stage2`, created from `main` at `69c2b4c` (the Stage-1 state is tagged `stage1-fallback`). This plan is the branch's first commit. **Never push** and never merge to `main` in this plan (both are gates for later).
- Baseline: `python -m pytest -q` → `308 passed` with a local pipeline configuration present (`305 passed, 3 skipped` without one). The Stage-1 golden test (`tests/test_stage1_golden.py`) must pass after every task.
- Schedule: subagents execute from 2026-10-15. The grillage (Tasks 1–8) must be finished by 2026-10-22, the student's exam day; B5 (Tasks 9–12) follows; due 2026-11-01. The student answers short gate questions only and reviews after the exam.
- **Gates in this plan:** Task 8 installs `PyNiteFEA==2.0.2` (approved 2026-10-07, see `docs/DECISIONS.md`). If the install or `import Pynite` fails, stop and report **BLOCKED** with the error output; do not substitute another package or version. Nothing else needs a student decision.
- Plan-writing probe: every code block below was run on a scratch clone of `69c2b4c` on 2026-10-15 (Python 3.13.3). The 68 new tests of Tasks 1–7 and 9–11 passed together with the whole Stage-1 suite (including the golden test), and `scripts/run_stage2.py` ran end to end in about 5.5 minutes without its PyNite part. PyNite was not installed then: the code in Task 8 follows the 2.0.2 source (`FEModel3D`, `Member3D`, `Node3D`) but has not been run. The code blocks of this file were then replayed onto a fresh clone of `69c2b4c` exactly as the steps say (appends, replacements, new files): `python -m pytest -q --ignore=tests/test_pynite_check.py` → `373 passed, 3 skipped` (no pipeline configuration on that clone). Numbers quoted from the probe are marked "probe"; the results that go into the log come only from the execution run.
- Where Stage-1 code is reused (not copied): `elements.beam_k`, `udl_equivalent`, `point_equivalent` (bending block and loads of the grid element); `girder1d.Girder1D` (one per girder line: staged sections `sec_el`/`sec_nd`, stage-A moments `M_A_el`/`M_A_nd`, `mcr`, `crack_at`, `eps_p0`, `EI_I`, `make_section`); `staging.design_prestress_staged` and `prestress_staged` (strand design and its check); `stiffness.ec2_secant_EI` (the same secant rule); `loads.truck_position_max_midspan`, `axle_offsets`, `impact_factor`, `truck_fraction`; `inverse.estimate_eta`, `eta_decompression`; `uq.monte_carlo_eta` (with the table as forward map); `claimbridge.governing_width`, `refused_below_resolution`; `figures.fig_mesh`. Stage-1 files are only appended to (`params.py`, `elements.py`, `nonlinear.py`, `figures.py`); no existing line changes.

## Conventions (Plan 1's conventions hold; Stage 2 adds these)

- **Units:** SI inside the code. Inputs whose name ends in `_ft` hold the source's feet and are converted once with `params.FT`.
- **Grillage axes:** x along the span (0 at the first support), z across the deck (girder G1 on z = 0, G5 on z = 4S), y up. Girder lines are indexed from 0: G1 = 0, G3 = 2, G5 = 4. Node (g, k) sits at (x_k, z_g) and has index g·(n_el + 1) + k.
- **Grid DOFs per node** (Logan Eq. 5.4.17 order): (w, θx, θz) → global (3i, 3i + 1, 3i + 2). w is the vertical deflection, up positive; θx and θz are rotations about the global x and z axes by the right-hand rule. For a girder member (along +x) θz = dw/dx, the beam's θ of Stage 1, and θx is its twist. For a cross member (along +z) the bending slope is dw/dz = −θx and θz is its twist.
- **Element:** local DOFs (w1, φ1x, φ1z, w2, φ2x, φ2z) with x̂ from node i to node j, ŷ = y, ẑ = x̂ × ŷ. The bending block (w, φz) is `beam_k`; the torsion block is GJ/L·[[1, −1], [−1, 1]]. **[T]** (Logan Eq. 5.4.18) per node is [[1, 0, 0], [0, C, S], [0, −S, C]] with C = (xj − xi)/L, S = (zj − zi)/L; local = T·global, k = Tᵀ k̂ T (Eq. 5.4.19), global consistent loads = Tᵀ f̂_eq.
- **Member end forces** are local: f̂ = k̂ T d − f̂_eq = [V1, T1, M1, V2, T2, M2] acting on the member. Sagging moment at x̂ from node i: M(x̂) = −M1 + V1·x̂ + q·x̂²/2 + Σ P·(x̂ − a), exactly the Stage-1 rule. **Torsion sign:** the internal torque is −T1 (Logan Eq. 5.4.9), constant along the member; a torque about +x is positive by the right-hand rule.
- **Loads off the girder lines:** shared by the lever rule between the two adjacent girder lines; on the overhang the whole load goes to the exterior girder with a torque m_x = (z_load − z_girder)·P for a downward load of magnitude P. Line torques become consistent nodal moments (half of each element's share to each end node).
- **Moment at a girder node:** cross members deliver moment into a girder at its nodes, so the girder moment jumps there; `line_node_moments` reports the mean of the two sides at an interior node.
- **Stage 2 keeps Stage 1's physics conventions:** `f_eq = +{f0}`; staged sections carry the slab's locked-in strain `eps_free`; secant EI = (M − M_ref)/(κ(M) − κ(M_ref)) with M_ref = the element's own stage-A moment; η scales today's strand strain; Monte Carlo σ = u/2, truncated at w > 0, seed `uq.DEFAULT_SEED`; deflection bands come from η percentiles in reversed order.
- **Stage A vs stage B in the grillage:** stage A (prestress, girder, wet slab over the girder's own slab width) acts on each bare girder and comes from statics; stage B (wearing surface, barriers, truck) acts on the grillage. The cracked girder Gk has η < 1, every other girder η = 1, all under the same secant rule.
- **Mesh:** `n_el = 20` girder elements along the span (Task 12 records the convergence study 10–80); cross members at every girder station; truck placements on a 0.05 m transverse step.
- **Commands** run from `fem_project/` (Git Bash). Tests: `python -m pytest ...`. Every commit ends with the trailer `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Repo hygiene:** no absolute local paths, no ClaimBridge internals (code, CLI names, internal paths) in committed text; the pipeline command lives only in the ignored `pipeline.local.json`. Cite the textbook as "Logan (5th ed.) §x, p. N", never a file path.

## File structure

| File | Status | Responsibility |
|---|---|---|
| `src/crackfem/params.py` | append | Stage-2 inputs (overhang, barrier width, wheel/lane geometry, ν, diaphragm), each with a source id |
| `src/crackfem/elements.py` | append | grid element: local k, [T], global k, consistent loads (Logan §5.4) |
| `src/crackfem/gridmodel.py` | create | grid model data: nodes (x, z), members, supports, member/nodal loads |
| `src/crackfem/gridsolve.py` | create | assembly, partition solve with a reusable Cholesky factor, member end forces, moments, torque |
| `src/crackfem/nonlinear.py` | append | `secant_iterate`: generic secant loop returning `SecantResult` (converged flag, history) |
| `src/crackfem/deck.py` | create | S01 deck layout, Hambly member properties, lever-rule load mapping, stage-B dead load, wheel loads, truck placements |
| `src/crackfem/grillage.py` | create | Stage-2 wiring: per-girder `Girder1D`, strand design, critical placements, truck envelope, nonlinear analysis, w_pred, indicator |
| `src/crackfem/pynite_check.py` | create | the same linear grid solved in PyNite |
| `src/crackfem/wtable.py` | create | tabulated w(η) with the plateau edge, monotone interpolation, error check |
| `src/crackfem/figures.py` | append | fig9–fig14 drawing functions |
| `scripts/run_stage2.py` | create | all Stage-2 results → `results/stage2.json` + figures |
| `tests/test_params_stage2.py`, `test_grid_element.py`, `test_gridsolve.py`, `test_secant_iterate.py`, `test_deck.py`, `test_grillage_linear.py`, `test_grillage.py`, `test_grillage_nonlinear.py`, `test_pynite_check.py`, `test_wtable.py`, `test_figures_stage2.py`, `test_run_stage2.py` | create | one file per module or check |
| `tests/fixtures/published.json` | modify | add `logan_grid_example` (Logan Example 5.5, printed values) |
| `requirements.txt` | modify | pin `PyNiteFEA==2.0.2` and its new transitive packages |
| `docs/sources.md` | modify | fill S22; add S23–S27 |
| `docs/DECISIONS.md`, `docs/PROGRESS_LOG.md`, `README.md` | modify | decisions, Entry 4, Stage-2 quickstart |
| `results/stage2.json`, `results/figures/fig9_*.png` … `fig14_*.png` | generated | Stage-2 results (Stage-1 outputs untouched) |

## Task groups (one subagent per group; review after each group)

| Group | Tasks | What exists afterwards | Review focus |
|---|---|---|---|
| B4-a | 1–4 | inputs and sources, grid element (sympy + Logan Ex. 5.5), grid solve, generic secant loop | element signs and [T]; Logan numbers from the fixture, not from the code |
| B4-b | 5–7 | deck model, linear grillage checks, `Grillage` (design, placements, envelope, nonlinear loop, w_pred, indicator) | statics of the load mapping; zero-tension design; loop > 1 pass |
| B4-c | 8 | PyNite installed, pinned and agreeing with the grid solve | **gate**: install; BLOCKED if it fails |
| B5 | 9–12 | table, figures, run script, `results/stage2.json`, Entry 4, decisions, README | numbers in the log trace to `stage2.json`; Stage-1 outputs unchanged |

---

## Task 1: Stage-2 inputs, sources S22–S27, design decisions

**Files:**
- Modify: `src/crackfem/params.py` (append fields to `Inputs` after `deflection_limit_ratio`; append entries to `SOURCES`)
- Modify: `docs/sources.md` (replace the S22 row; add rows S23–S27 after it)
- Modify: `docs/DECISIONS.md` (append a section)
- Test: `tests/test_params_stage2.py`

- [ ] **Step 1: Write the failing test** — `tests/test_params_stage2.py`

```python
"""Stage-2 inputs: values and source rows (S22-S27)."""
import re
from pathlib import Path

import pytest

from crackfem.params import FT, SOURCES, Inputs

SOURCES_MD = Path(__file__).resolve().parents[1] / "docs" / "sources.md"
STAGE2 = {"overhang_width_ft": (3.0, "S25"), "barrier_width_ft": (1.0, "S25"),
          "wheel_spacing_ft": (6.0, "S26"), "truck_clearance_ft": (10.0, "S26"),
          "wheel_curb_clearance_ft": (2.0, "S26"), "lane_width_ft": (12.0, "S26"),
          "poisson_concrete": (0.2, "S24"), "diaphragm_depth": (1.10, "S27"),
          "diaphragm_thickness": (0.25, "S27")}


@pytest.mark.parametrize("name", sorted(STAGE2))
def test_stage2_input_value_and_source(name):
    value, sid = STAGE2[name]
    assert getattr(Inputs(), name) == value
    assert SOURCES[name] == sid


def test_new_source_rows_exist_and_s22_is_filled():
    text = SOURCES_MD.read_text(encoding="utf-8")
    ids = re.findall(r"^\|\s*(S\d\d)\s*\|", text, flags=re.MULTILINE)
    assert {"S22", "S23", "S24", "S25", "S26", "S27"} <= set(ids)
    s22 = next(line for line in text.splitlines() if line.startswith("| S22 |"))
    assert "Example 5.5" in s22 and "local PDF" not in s22


def test_truck_fits_its_clearance_width():
    inp = Inputs()
    assert 2 * inp.wheel_curb_clearance_ft + inp.wheel_spacing_ft == pytest.approx(inp.truck_clearance_ft)
    assert inp.overhang_width_ft * FT > inp.barrier_width_ft * FT
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_params_stage2.py -q`
Expected: FAIL — `AttributeError: 'Inputs' object has no attribute 'barrier_width_ft'` (and the other new names), and the S22–S27 row check fails.

- [ ] **Step 3: Append the Stage-2 inputs to `src/crackfem/params.py`**

In class `Inputs`, directly after the line `    deflection_limit_ratio: float = 800.0  # span / limit`, add:

```python
    # grillage (Stage 2)
    overhang_width_ft: float = 3.0         # ft, exterior girder centreline to deck edge
    barrier_width_ft: float = 1.0          # ft, deck edge to the curb face of the barrier
    wheel_spacing_ft: float = 6.0          # ft, transverse distance between the two wheel lines
    truck_clearance_ft: float = 10.0       # ft, clearance width of one truck (wheels 2 ft inside)
    wheel_curb_clearance_ft: float = 2.0   # ft, smallest distance from a wheel line to the curb face
    lane_width_ft: float = 12.0            # ft, design lane width (number of lanes = floor(w / 12))
    poisson_concrete: float = 0.2          # uncracked concrete, dimensionless
    diaphragm_depth: float = 1.10          # m, diaphragm web depth below the slab
    diaphragm_thickness: float = 0.25      # m, diaphragm web thickness
```

Replace the whole `SOURCES = {...}` assignment with (the Stage-1 entries unchanged, five Stage-2 lines added at the end):

```python
SOURCES = {
    "span": "S01",
    "girder_spacing": "S02", "slab_thickness": "S02", "wearing_thickness": "S02",
    "n_girders": "S02", "n_barriers": "S02", "barrier_weight": "S02",
    "girder_outline_in": "S03",
    "strand_diameter": "S04", "strand_area": "S04", "fpu": "S04",
    "Ep": "S05",
    "fpj_ratio": "S06", "lump_sum_loss": "S06",
    "strand_centroid": "S15", "strand_cover": "S15",
    "fck_girder": "S19", "fck_slab": "S19",
    "unit_weight_concrete": "S11", "unit_weight_wearing": "S11",
    "k1": "S08", "k2": "S08", "k3": "S08", "k4": "S08", "kt": "S08",
    "bond_xi": "S09",
    "beta": "S10",
    "axle_loads": "S12", "axle_spacings": "S12",
    "impact_numerator_ft": "S13", "impact_denominator_ft": "S13", "impact_cap": "S13",
    "df_divisor_ft": "S14",
    "live_load_factor": "S12",
    "deflection_limit_ratio": "S18",
    "overhang_width_ft": "S25", "barrier_width_ft": "S25",
    "wheel_spacing_ft": "S26", "truck_clearance_ft": "S26", "wheel_curb_clearance_ft": "S26",
    "lane_width_ft": "S26",
    "poisson_concrete": "S24",
    "diaphragm_depth": "S27", "diaphragm_thickness": "S27",
}
```

- [ ] **Step 4: Replace the S22 row of `docs/sources.md` and add S23–S27**

Replace the whole line `| S22 | Logan §5.4 grid worked example (Plan 2) | fixture | local PDF, Plan 2 | | |` with these six lines (one table row each):

```markdown
| S22 | Grid element (local k, [T], k = TᵀkT), torsion constant of an open section, and the worked grid example for the tests | fixture `logan_grid_example` (Example 5.5: E 30,000 ksi, G 12,000 ksi, I 400 in⁴, J 110 in⁴, 100 kip at node 1); Type IV J = Σ b t³/3 over the bottom flange 26 × 8, top flange 20 × 8 and web 38 × 8 in (haunches left out) = 14,336 in⁴ = 0.005967 m⁴ | Logan, *A First Course in the Finite Element Method*, 5th ed., §5.4 Grid Equations: Eq. 5.4.15 (J ≈ Σ b t³/3, open sections) p. 241; Table 5-1 item 4 (wide-flange beam with unequal flanges) p. 242; Eq. 5.4.17 (local k, DOF order deflection, torsional rotation, bending rotation) p. 243; Eq. 5.4.18 ([T], C = (xj − xi)/L, S = (zj − zi)/L) p. 243; Eq. 5.4.19 p. 244; Example 5.5 pp. 244–250: printed K (Eq. 5.4.27), d1y = −2.83 in, φ1x = 0.0295 rad, φ1z = −0.0169 rad (Eq. 5.4.29), local forces Eqs. 5.4.32, 5.4.34, 5.4.36. | [primary] | 2026-10-15 |
| S23 | Grillage member properties (Hambly's rules) | slab strip of width b: I = b t³/12, J = b t³/6; girder: J = J_TypeIV + b t³/6 with b its slab width (sum of girder and slab, St Venant); exterior girder slab width S/2 + overhang | Yang, Miao & Ghosn, *Redundancy Analysis of Composite Spread Box Girder Superstructures under Vertical Loads*, NCHRP Report 776, Appendix D.1, §2.1, p. D.1-3, Eq. (2.1): "J = t³/6 per unit length of slab"; as proposed by Hambly (1991), half the thin-rectangle value t³/3 to account for slab continuity — <https://onlinepubs.trb.org/onlinepubs/nchrp/nchrp_rpt_776AppendixD1.pdf>. LUSAS, *Simple Grillage* worked example (version 21.0), pp. 4 and 6: grillage formulae "taken from Hambly's Bridge Deck Behaviour", summarised in CSA S6-14 Annex A5.1; default torsion setting "Sum slab and girder", St Venant — <https://www.lusas.com/user_area/documentation/V21_0/worked_examples/Simple%20Grillage.pdf>. I = b t³/12 is the second moment of a rectangle. Hambly (*Bridge Deck Behaviour*) itself was not opened. Exterior slab width checked against S17: 2.01 m, each side within 6 t + web/2 = 1.30 m. | [secondary] | 2026-10-15 |
| S24 | Poisson's ratio of uncracked concrete (shear modulus G = E / 2(1 + ν)) | 0.2 | *Table of concrete design properties (fcd, fctm, Ecm, fctd) – Eurocode 2*, eurocodeapplied.com: "ν = 0.2 for uncracked concrete", EN 1992-1-1 §3.1.3(4); G = E / [2 (1 + ν)] in the same table — <https://cdn.eurocodeapplied.com/design/en1992/concrete-design-properties>. EN 1992-1-1 itself was not opened. | [secondary] | 2026-10-15 |
| S25 | Deck overhang and barrier width | overhang 3.0 ft (0.914 m) from the exterior girder centreline; barrier (curb face to deck edge) 1.0 ft (0.305 m) | TTI 0-4751-1 Vol. 2 (see S02), §A.2.5.2.2.1, p. A.2-11: de = 3.0 − 0.33 − 1.0 = 1.67 ft, i.e. overhang 3.0 ft from the centre of the exterior girder and a 1.0 ft barrier width (T501, the barrier of S02). Same example as S02; not the real CPK72 deck. | [stand-in] | 2026-10-15 |
| S26 | HS truck across the deck and design lanes | wheel lines 6 ft apart, 2 ft inside a 10 ft clearance width placed against the curb at the closest; lanes 12 ft, number of design lanes = floor(roadway / 12 ft) | Ajeel, *Structural Analysis of Concrete Bridges – Loads on Bridge Structures* (see S12), Fig. 2.5 (H20-44 and HS20-44): "2'-0" 6'-0" 2'-0"" across, "Clearance and Lane Width 10'-0"", curb shown. TTI 0-4751-1 Vol. 2, p. A.2-11: "Number of design lanes = Integer part of the ratio w/12", w = clear roadway width between the curbs [LRFD Art. 3.6.1.1.1]. The AASHTO Standard Specifications' lane article itself was not opened. | [secondary] | 2026-10-15 |
| S27 | Diaphragms: end and midspan, depth 1.10 m below the slab, thickness 0.25 m | as listed; modelled as the web alone (not composite with the slab), added to the slab strip of its station | No source: a cast-in-place diaphragm spanning most of the 1.37 m Type IV depth below the slab. Placement (both ends and midspan) is the spec §5.1 deck layout (S02). Not the real CPK72 deck. | [stand-in] | 2026-10-15 |
```

- [ ] **Step 5: Append the design decisions to `docs/DECISIONS.md`**

```markdown
## 2026-10-15 — Plan 2 design (Stage-2 grillage, B4–B5; controller-approved list, plan-writing probe)

- **Grid element = Logan §5.4** (Eqs. 5.4.16–5.4.19): 3 DOF per node in Logan's order (w, θx, θz), bending block = the Stage-1 `beam_k` on (w, θz), torsion GJ/L, [T] per node [[1, 0, 0], [0, C, S], [0, −S, C]]; checked against a sympy derivation and Logan (5th ed.) Example 5.5 (pp. 244–250, printed answers in `tests/fixtures/published.json`). Rejected: a 3-D frame element (the in-plane DOFs of a grid carry nothing).
- **Deck model:** 5 girder lines at 2.2 m with the same longitudinal mesh as the girders; slab strips at every girder station with tributary width (Hambly: I = b t³/12, J = b t³/6, S23); end and midspan diaphragms; overhang. Rejected: strips at a coarser spacing than the girder nodes (loads and moments would need interpolation between stations).
- **Mesh n_el = 20:** probe convergence study (strand counts held): w_pred at η = 0.36 within 0.03 % of N = 80 for G1 and G3; damaged service deflection of G1 within 0.6 % (the moving-truck envelope that sets the pre-cracked zone is resolved on the element grid). N = 40 costs 4× the build time for no visible change. The study is rerun in Task 12.
- **Diaphragm = its web alone** (1.10 × 0.25 m, slab concrete, S27) added to the slab strip of its station, not a T section with the strip as flange. Reason: the probe showed the T-section diaphragm's stiffness depends on the mesh (the flange is the strip width L/n_el), and G1's w_pred drifted by 3 % per mesh doubling with it; with the web alone the result is mesh-independent. A non-composite diaphragm is a lower bound on its stiffness (stated).
- **Overhang 3.0 ft and barrier width 1.0 ft** from the TxDOT/TTI Type IV example (S25, the source of S02), instead of a free 1.0 m: the same example family as the other deck stand-ins. With the 2 ft curb rule the outer wheel line then lies exactly over G1 at the curb position.
- **Member properties:** girder longitudinal EI from the Stage-1 section models of a composite girder with its own slab width (exterior S/2 + overhang = 2.01 m, interior S = 2.2 m); girder GJ = G_girder × J_TypeIV + G_slab × b t³/6 (girder plus slab, St Venant; J_TypeIV = Σ b t³/3 over bottom flange, top flange and web, haunches left out, Logan Eq. 5.4.15 and Table 5-1, = 14,336 in⁴); G = Ecm/2(1 + ν) with ν = 0.2 (EC2 3.1.3(4), S24). Cross members and torsion stay linear (uncracked).
- **Stage A per girder by statics** exactly as Stage 1, each `Girder1D` built with `girder_spacing` = its slab width, so the exterior girder carries its own girder and slab tributary. The wet-slab overhang's torque on the bare exterior girder is not modelled (stage A is vertical statics).
- **Stage-B dead load as girder line loads by the lever rule**, cut at the girder lines (exact force and moment), with the overhang part and the barriers as a line torque on the exterior girders turned into consistent nodal moments. Rejected: nodal loads at the stations (statically equal, but the moment inside each girder element would lose its parabolic part).
- **Truck: wheel loads, no distribution factor.** HS20-44 × 1.30 × (1 + I), each axle split into two wheel lines 6 ft apart, wheels 2 ft inside a 10 ft clearance width that stays between the curbs (S26). Wheel lines between girders go to the two adjacent girder lines by the lever rule, then act as point loads on the girder members (consistent loads). Longitudinal position = the Stage-1 placement for maximum midspan moment.
- **Deviation from "the truck" (singular): one truck and two trucks side by side are scanned** (two design lanes: roadway 10.02 m, floor(w / 12 ft) = 2; clearance widths may not overlap; every placement also mirrored). Reason: Stage 1's S/5.5 is the AASHTO factor for two or more lanes loaded, so a one-truck grillage would not be comparable; in the probe two trucks govern G1. The critical placement of each girder (largest midspan stage-B truck moment on the linear intact grillage, by superposition of one solve per girder line) is held fixed for every η: a crack does not move the truck.
- **Nonlinear grillage:** every girder element gets the Stage-1 rule (staged section, EC2 secant EI from M_ref = its own stage-A moment); the cracked girder Gk has η < 1 and the others η = 1, so a neighbour that cracks under the shed load is handled by the same rule (the probe found G2 cracking when G1 is at η = 0.36). Sticky pre-cracked flag inside an analysis; stiffness-residual check; under-relaxation available (`relax`); iteration cap 100; a non-converged analysis is returned with its history (`SecantResult.converged = False`) and makes w_pred NaN, which the inverse refuses. Sensitivity: GJ of cracked girder elements × 0.5 (one run per cracked girder at its reference η).
- **Initial pre-cracked zone from a linear intact envelope:** stage A + stage-B dead + the moving truck at the girder's critical placement (each axle at every node and element midpoint, both directions), compared with Mcr(η). Rejected: a nonlinear moving-truck envelope per η (hundreds of nonlinear solves per w_pred). The envelope over-estimates the cracked girder's truck moment (it ignores its own shedding), so it errs toward more pre-cracking.
- **Strand design per girder:** the Stage-1 zero-tension rule under each girder's own grillage service moment at η = 1 (stage A + stage-B dead + its critical truck), iterated until the counts stop changing; probe (31, 29, 27, 29, 31) in two rounds. Rejected: 30 strands everywhere (G1 would have soffit tension at service, the interior girders a different margin than the rule gives, and η would not mean the same thing for G1 and G3).
- **Cracked girders:** G1 (index 0, exterior, default) and G3 (index 2, interior), as in the spec §5.4.
- **Inverse and Monte Carlo through a table:** w_pred on a 99-point η grid (0.02–1) plus the plateau edge (bisection on the exact map), linear interpolation; `uq.monte_carlo_eta` (2,000 samples, σ = u/2, same seed) runs on the table unchanged; p05/p50/p95 are confirmed by the exact `estimate_eta` at the widths the table maps them to; the point estimate is exact. Rejected: Monte Carlo on the exact map (about 25,000 nonlinear solves per case). The data (`data/claimbridge/*_spike.json`) and the inverse/uq code are reused unchanged.
- **Reported moments are stage-B (grillage-carried) midspan moments**; stage A comes from statics and is the same intact or damaged. The moment at a girder node is the mean of the two sides (cross members deliver a moment jump at the node).
- **Outputs:** `results/stage2.json` and fig9–fig14; Stage-1 outputs are read, never written; the golden test stays. **Deviation from spec §8 (B5): no tornado on the grillage** — the Plan-2 scope set on 2026-10-15 lists the outputs without it; the Stage-1 tornado (as-built strands) stands.
- **PyNiteFEA pinned at 2.0.2**, the newest release that does not require numpy ≥ 2.4 (2.1.0–3.2.0 declare `numpy>=2.4.0`, which would move the pinned numpy 2.3.5 and possibly the Stage-1 numbers); installed with numpy, scipy and matplotlib held at their pins. New transitive package: PrettyTable.
```

- [ ] **Step 6: Run the tests**

Run: `python -m pytest tests/test_params_stage2.py tests/test_params.py tests/test_stage1_golden.py -q`
Expected: `19 passed` (11 + 3 + 5).

- [ ] **Step 7: Commit**

```bash
git add src/crackfem/params.py docs/sources.md docs/DECISIONS.md tests/test_params_stage2.py
git commit -m "feat: Stage-2 inputs with sources S22-S27, Plan 2 decisions" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 2: Grid element (Logan §5.4) + the Logan Example 5.5 fixture

**Files:**
- Modify: `src/crackfem/elements.py` (append; `beam_k`, `udl_equivalent`, `point_equivalent` reused)
- Modify: `tests/fixtures/published.json` (add the key `logan_grid_example`)
- Test: `tests/test_grid_element.py`

The worked example is Logan (5th ed.) §5.4, Example 5.5, pp. 244–250: three grid elements meeting at node 1, fixed at nodes 2, 3 and 4, a downward 100 kip load at node 1; E = 30,000 ksi, G = 12,000 ksi, I = 400 in⁴, J = 110 in⁴ for all elements; element lengths 22.36 ft, 22.36 ft and 10 ft. Printed answers: the reduced stiffness (Eq. 5.4.27), d1y = −2.83 in, φ1x = 0.0295 rad, φ1z = −0.0169 rad (Eq. 5.4.29), and the local element forces (Eqs. 5.4.32, 5.4.34, 5.4.36). The fixture keeps the printed US units; the test converts them. Logan works longhand to three significant figures (e.g. element 1 m̂2z printed −2660 k-in; exact −2652), so the tolerance is 0.5 %.

- [ ] **Step 1: Add the fixture key**

In `tests/fixtures/published.json`, add this key as the last member of the top-level object (put a comma after the closing brace of `pc_cracked_example`):

```json
  "logan_grid_example": {
    "source": "S22",
    "reference": "Logan (5th ed.) §5.4, Example 5.5, pp. 244-250",
    "units": "kip, in, ksi as printed; node coordinates in ft as in Figure 5-18",
    "E_ksi": 30000.0,
    "G_ksi": 12000.0,
    "I_in4": 400.0,
    "J_in4": 110.0,
    "nodes_ft": {
      "1": [
        20.0,
        10.0
      ],
      "2": [
        0.0,
        20.0
      ],
      "3": [
        0.0,
        0.0
      ],
      "4": [
        20.0,
        0.0
      ]
    },
    "elements": [
      [
        1,
        2
      ],
      [
        1,
        3
      ],
      [
        1,
        4
      ]
    ],
    "fixed_nodes": [
      2,
      3,
      4
    ],
    "load": {
      "node": 1,
      "F_y_kip": -100.0
    },
    "K_free_kip_in": [
      [
        98.2,
        5000.0,
        -1790.0
      ],
      [
        5000.0,
        479000.0,
        0.0
      ],
      [
        -1790.0,
        0.0,
        299000.0
      ]
    ],
    "displacements": {
      "d1y_in": -2.83,
      "phi1x_rad": 0.0295,
      "phi1z_rad": -0.0169
    },
    "local_forces_kip_kin": {
      "1": [
        -19.2,
        -167.0,
        -2480.0,
        19.2,
        167.0,
        -2660.0
      ],
      "2": [
        7.23,
        -92.5,
        2240.0,
        -7.23,
        92.5,
        -295.0
      ],
      "3": [
        -88.1,
        186.0,
        -2340.0,
        88.1,
        -186.0,
        -8240.0
      ]
    },
    "note": "Printed answers: K Eq. 5.4.27, displacements Eq. 5.4.29, local forces Eqs. 5.4.32, 5.4.34, 5.4.36 in the order [f1y, m1x, m1z, f2y, m2x, m2z] (element 2 ends at node 3, element 3 at node 4). Longhand to three significant figures, so the test tolerance is 0.5 %."
  }
```

- [ ] **Step 2: Write the failing test** — `tests/test_grid_element.py`

```python
"""Grid element (Logan 5th ed. §5.4): symbolic derivation and the printed Example 5.5."""
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from crackfem.elements import (beam_k, grid_k_global, grid_k_local, grid_point_local, grid_T,
                               grid_udl_local)
from crackfem.params import FT, INCH, KIP, KSI

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))
x, L, EI, GJ, q, P, a = sp.symbols("x L EI GJ q P a", positive=True)
xi = x / L
NW = sp.Matrix([1 - 3 * xi**2 + 2 * xi**3, L * (xi - 2 * xi**2 + xi**3), 3 * xi**2 - 2 * xi**3, L * (-xi**2 + xi**3)])
NT = sp.Matrix([1 - xi, xi])                         # linear twist, Logan Eq. 5.4.4


def _symbolic_k():
    kb = sp.integrate(EI * NW.diff(x, 2) * NW.diff(x, 2).T, (x, 0, L))
    kt = sp.integrate(GJ * NT.diff(x) * NT.diff(x).T, (x, 0, L))
    k = sp.zeros(6, 6)
    for r, R in enumerate((0, 2, 3, 5)):
        for c, Cc in enumerate((0, 2, 3, 5)):
            k[R, Cc] = kb[r, c]
    for r, R in enumerate((1, 4)):
        for c, Cc in enumerate((1, 4)):
            k[R, Cc] = kt[r, c]
    return sp.lambdify((EI, GJ, L), k)


def test_local_stiffness_equals_bending_plus_torsion_energy():
    k_sym = _symbolic_k()
    for EIv, GJv, Lv in ((2.0e9, 3.0e8, 1.3), (5.5, 0.7, 0.4)):
        np.testing.assert_allclose(grid_k_local(EIv, GJv, Lv), np.array(k_sym(EIv, GJv, Lv), float), rtol=1e-12)


def test_bending_block_is_the_beam_element():
    k = grid_k_local(3.0e9, 1.0e8, 2.0)
    idx = (0, 2, 3, 5)
    np.testing.assert_allclose(k[np.ix_(idx, idx)], beam_k(3.0e9, 2.0), rtol=1e-15)
    assert k[1, 1] == pytest.approx(1.0e8 / 2.0) and k[1, 4] == pytest.approx(-1.0e8 / 2.0)


def test_consistent_loads_are_the_work_integrals():
    fq = sp.lambdify((q, L), sp.integrate(q * NW, (x, 0, L)))
    fP = sp.lambdify((P, a, L), P * NW.subs(x, a))
    np.testing.assert_allclose(grid_udl_local(-12.5, 1.7)[[0, 2, 3, 5]], np.array(fq(-12.5, 1.7), float).ravel(), rtol=1e-12)
    np.testing.assert_allclose(grid_point_local(-3.0, 0.45, 1.2)[[0, 2, 3, 5]],
                               np.array(fP(-3.0, 0.45, 1.2), float).ravel(), rtol=1e-12)
    assert grid_udl_local(-12.5, 1.7)[1] == 0.0 and grid_point_local(-3.0, 0.45, 1.2)[4] == 0.0


@pytest.mark.parametrize("theta", [0.0, 0.3, 1.2, np.pi / 2, 2.5])
def test_transformation_is_orthogonal_and_keeps_the_deflection(theta):
    T = grid_T(np.cos(theta), np.sin(theta))
    np.testing.assert_allclose(T @ T.T, np.eye(6), atol=1e-15)
    assert T[0, 0] == 1.0 and T[3, 3] == 1.0


def test_global_stiffness_is_symmetric_with_rigid_body_modes():
    k = grid_k_global(2.0e9, 3.0e8, 2.5, 0.6, 0.8)
    np.testing.assert_allclose(k, k.T, rtol=1e-12, atol=1e-3)
    rigid = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])       # vertical translation
    np.testing.assert_allclose(k @ rigid, 0.0, atol=1e-3)


def _logan_example():
    ex = FIX["logan_grid_example"]
    E, G = ex["E_ksi"] * KSI, ex["G_ksi"] * KSI
    I, J = ex["I_in4"] * INCH**4, ex["J_in4"] * INCH**4
    xz = {int(n): np.array(c) * FT for n, c in ex["nodes_ft"].items()}
    return ex, E * I, G * J, xz


def test_logan_example_5_5_stiffness_and_displacements():
    ex, EIv, GJv, xz = _logan_example()
    K = np.zeros((3, 3))
    for i, j in ex["elements"]:
        d = xz[j] - xz[i]
        Lm = float(np.hypot(*d))
        K += grid_k_global(EIv, GJv, Lm, d[0] / Lm, d[1] / Lm)[:3, :3]
    # K in kip/in, kip, kip-in: rows/cols (w, phi_x, phi_z); w in in, phi in rad
    to_us = np.array([[KIP / INCH, KIP, KIP], [KIP, KIP * INCH, KIP * INCH], [KIP, KIP * INCH, KIP * INCH]])
    np.testing.assert_allclose(K / to_us, ex["K_free_kip_in"], rtol=5e-3, atol=1e-6)
    d = np.linalg.solve(K, [ex["load"]["F_y_kip"] * KIP, 0.0, 0.0])
    p = ex["displacements"]
    np.testing.assert_allclose([d[0] / INCH, d[1], d[2]], [p["d1y_in"], p["phi1x_rad"], p["phi1z_rad"]], rtol=5e-3)


def test_logan_example_5_5_local_element_forces():
    ex, EIv, GJv, xz = _logan_example()
    K = np.zeros((3, 3))
    geo = []
    for i, j in ex["elements"]:
        dd = xz[j] - xz[i]
        Lm = float(np.hypot(*dd))
        geo.append((Lm, dd[0] / Lm, dd[1] / Lm))
        K += grid_k_global(EIv, GJv, *geo[-1])[:3, :3]
    d1 = np.linalg.solve(K, [ex["load"]["F_y_kip"] * KIP, 0.0, 0.0])
    to_us = np.array([KIP, KIP * INCH, KIP * INCH, KIP, KIP * INCH, KIP * INCH])
    for e, (Lm, C, S) in enumerate(geo, start=1):
        f = grid_k_local(EIv, GJv, Lm) @ grid_T(C, S) @ np.concatenate([d1, np.zeros(3)])
        np.testing.assert_allclose(f / to_us, ex["local_forces_kip_kin"][str(e)], rtol=5e-3)
```

- [ ] **Step 3: Run it to verify it fails**

Run: `python -m pytest tests/test_grid_element.py -q`
Expected: FAIL — collection error `ImportError: cannot import name 'grid_k_global' from 'crackfem.elements'`.

- [ ] **Step 4: Append the grid element to `src/crackfem/elements.py`**

Add at the end of the file (after `point_equivalent`):

```python
# --- Grid element (Logan 5th ed. §5.4, Eqs. 5.4.16-5.4.19) -------------------------------
# Local DOFs per element: (w1, phi1x, phi1z, w2, phi2x, phi2z): vertical deflection (y up +),
# twist about the member axis x^, bending rotation about z^ = x^ cross y^ (the beam's theta).
# The bending block is beam_k on (w, phi_z); the torsion block is GJ/L [[1, -1], [-1, 1]].
GRID_BENDING = (0, 2, 3, 5)
GRID_TORSION = (1, 4)


def grid_k_local(EI: float, GJ: float, L: float) -> np.ndarray:
    k = np.zeros((6, 6))
    k[np.ix_(GRID_BENDING, GRID_BENDING)] = beam_k(EI, L)
    t = GJ / L
    k[np.ix_(GRID_TORSION, GRID_TORSION)] = t * np.array([[1.0, -1.0], [-1.0, 1.0]])
    return k


def grid_T(C: float, S: float) -> np.ndarray:
    """Logan Eq. 5.4.18: local = T @ global per node, with C = (xj - xi) / L, S = (zj - zi) / L."""
    r = np.array([[1.0, 0.0, 0.0], [0.0, C, S], [0.0, -S, C]])
    T = np.zeros((6, 6))
    T[:3, :3] = r
    T[3:, 3:] = r
    return T


def grid_k_global(EI: float, GJ: float, L: float, C: float, S: float) -> np.ndarray:
    """Logan Eq. 5.4.19: k = T^T k^ T."""
    T = grid_T(C, S)
    return T.T @ grid_k_local(EI, GJ, L) @ T


def grid_udl_local(q: float, L: float) -> np.ndarray:
    """Consistent local loads of a uniform vertical load q (up +): the beam's, no torque."""
    f = np.zeros(6)
    f[list(GRID_BENDING)] = udl_equivalent(q, L)
    return f


def grid_point_local(P: float, a: float, L: float) -> np.ndarray:
    """Consistent local loads of a vertical point load P (up +) at distance a from node 1."""
    f = np.zeros(6)
    f[list(GRID_BENDING)] = point_equivalent(P, a, L)
    return f
```

- [ ] **Step 5: Run the tests**

Run: `python -m pytest tests/test_grid_element.py tests/test_elements.py -q`
Expected: `15 passed` (11 + 4).

- [ ] **Step 6: Commit**

```bash
git add src/crackfem/elements.py tests/fixtures/published.json tests/test_grid_element.py
git commit -m "feat: grid element (Logan 5.4) checked against sympy and Logan Example 5.5" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 3: Grid model and partition solve

**Files:**
- Create: `src/crackfem/gridmodel.py`, `src/crackfem/gridsolve.py`
- Test: `tests/test_gridsolve.py`

The course method (Ch3 slides; Logan §5.4): assemble k = Tᵀ k̂ T member by member, partition K_ff D_f = F_f, reactions R = K D − F at the restrained DOFs, member end forces f̂ = k̂ T d − f̂_eq. `GridStiffness` keeps the Cholesky factor of K_ff so that many load cases on one stiffness (truck placements, the moving-truck envelope) cost one back-substitution each.

- [ ] **Step 1: Write the failing test** — `tests/test_gridsolve.py`

```python
"""Grid assembly and partition solve: Logan Example 5.5 through the whole path, and the
reductions a grid must obey (a single line is a beam, a cantilever twists by TL/GJ)."""
import json
from pathlib import Path

import numpy as np
import pytest

from crackfem.gridmodel import GridModel, with_loads
from crackfem.gridsolve import line_node_moments, midpoint_moments, moment_at, solve, stiffness, torque
from crackfem.model import simply_supported
from crackfem.model import with_loads as beam_loads
from crackfem.params import FT, INCH, KIP, KSI
from crackfem.solve import node_moments
from crackfem.solve import solve as beam_solve

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))
L, EI, GJ, W = 24.5, 7.5e9, 1.2e9, 28e3


def test_logan_example_5_5_through_assembly_and_partition():
    ex = FIX["logan_grid_example"]
    names = sorted(ex["nodes_ft"], key=int)                       # "1".."4" -> indices 0..3
    xz = np.array([ex["nodes_ft"][n] for n in names]) * FT
    members = np.array(ex["elements"]) - 1
    fixed = tuple(3 * (n - 1) + d for n in ex["fixed_nodes"] for d in range(3))
    model = with_loads(GridModel(xz, members, fixed), nodal_loads=[(0, ex["load"]["F_y_kip"] * KIP)])
    E, G = ex["E_ksi"] * KSI, ex["G_ksi"] * KSI
    sol = solve(model, E * ex["I_in4"] * INCH**4, G * ex["J_in4"] * INCH**4)
    p = ex["displacements"]
    np.testing.assert_allclose([sol.D[0] / INCH, sol.D[1], sol.D[2]], [p["d1y_in"], p["phi1x_rad"], p["phi1z_rad"]], rtol=5e-3)
    to_us = np.array([KIP, KIP * INCH, KIP * INCH, KIP, KIP * INCH, KIP * INCH])
    for e in range(3):
        np.testing.assert_allclose(sol.end_forces[e] / to_us, ex["local_forces_kip_kin"][str(e + 1)], rtol=5e-3)
    assert sol.R[0::3].sum() == pytest.approx(-ex["load"]["F_y_kip"] * KIP, rel=1e-12)


def _line(n, along_z=False):
    s = np.linspace(0.0, L, n + 1)
    xz = np.stack([np.zeros_like(s), s], axis=1) if along_z else np.stack([s, np.zeros_like(s)], axis=1)
    members = np.array([(k, k + 1) for k in range(n)])
    twist = 2 if along_z else 1                                  # the DOF about the line's own axis
    return GridModel(xz, members, (0, twist, 3 * n))            # w at both ends, twist at one end


@pytest.mark.parametrize("along_z", [False, True])
def test_a_single_line_is_the_beam(along_z):
    n, pts = 8, [(5.0, -80e3), (17.0, -40e3)]
    beam = beam_loads(simply_supported(L, n), q=-W, point_loads=pts)
    ref = beam_solve(beam, EI)
    grid = _line(n, along_z)
    gl = [(int(xp // (L / n)), xp - (L / n) * int(xp // (L / n)), P) for xp, P in pts]
    grid = with_loads(grid, q=-W, point_loads=gl)
    sol = solve(grid, EI, GJ)
    np.testing.assert_allclose(sol.D[0::3], ref.D[0::2], rtol=1e-10, atol=1e-14)
    # bending rotation: theta_z for a line along x; along z the slope dw/dz is -theta_x
    rot = -sol.D[1::3] if along_z else sol.D[2::3]
    np.testing.assert_allclose(rot, ref.D[1::2], rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(line_node_moments(grid, sol, range(n)), node_moments(beam, ref), rtol=1e-9, atol=1e-3)
    np.testing.assert_allclose(midpoint_moments(grid, sol, range(n)),
                               [moment_at(grid, sol, m, L / n / 2) for m in range(n)], rtol=1e-12)


def test_cantilever_twists_by_TL_over_GJ():
    n, T = 4, 15e3
    s = np.linspace(0.0, 3.0, n + 1)
    model = GridModel(np.stack([s, np.zeros_like(s)], axis=1), np.array([(k, k + 1) for k in range(n)]), (0, 1, 2))
    sol = solve(with_loads(model, nodal_loads=[(3 * n + 1, T)]), EI, GJ)
    assert sol.D[3 * n + 1] == pytest.approx(T * 3.0 / GJ, rel=1e-12)
    assert all(torque(sol, m) == pytest.approx(T, rel=1e-12) for m in range(n))
    assert np.allclose(sol.D[0::3], 0.0)                                  # pure twist, no bending


def test_two_by_two_grid_is_symmetric_and_in_equilibrium():
    xz = np.array([(0, 0), (6, 0), (0, 4), (6, 4)], dtype=float)
    members = np.array([(0, 1), (2, 3), (0, 2), (1, 3)])
    fixed = (0, 3, 6)                                                    # three corners vertically
    model = with_loads(GridModel(xz, members, fixed + (1, 2)), q=[-5e3, -2e3, 0.0, -1e3],
                       point_loads=[(0, 2.0, -30e3)], nodal_loads=[(9, -10e3), (10, 4e3)])
    st = stiffness(model, [EI, EI, 0.5 * EI, 0.5 * EI], GJ)
    np.testing.assert_allclose(st.K, st.K.T, rtol=1e-12, atol=1e-3)
    sol = solve(model, stiff=st)
    total = 5e3 * 6 + 2e3 * 6 + 1e3 * 4 + 30e3 + 10e3
    assert sol.R[0::3].sum() == pytest.approx(total, rel=1e-12)


def test_reused_factorisation_gives_the_same_answer():
    m = with_loads(_line(6), q=-W)
    st = stiffness(m, EI, GJ)
    a, b = solve(m, EI, GJ), solve(m, stiff=st)
    np.testing.assert_array_equal(a.D, b.D)


def test_invalid_input_is_refused():
    with pytest.raises(ValueError):
        GridModel(np.zeros((2, 2)), np.array([(0, 1)]), (0,))            # zero-length member
    with pytest.raises(ValueError):
        solve(with_loads(_line(4), q=-W), -1.0, GJ)
    with pytest.raises(ValueError):
        with_loads(_line(4), point_loads=[(0, 99.0, -1.0)])
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_gridsolve.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'crackfem.gridmodel'`.

- [ ] **Step 3: Write `src/crackfem/gridmodel.py`**

```python
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
    """Copy of the model with these loads; replaces ALL existing loads. q: scalar or per member."""
    q = np.broadcast_to(np.asarray(q, dtype=float), (model.n_members,)).copy()
    L = model.lengths
    for m, a, _ in point_loads:
        if not 0 <= m < model.n_members or not 0.0 <= a <= L[m]:
            raise ValueError(f"point load at a = {a} is not on member {m}")
    for dof, _ in nodal_loads:
        if not 0 <= dof < model.n_dof:
            raise ValueError(f"nodal load DOF {dof} is out of range")
    return replace(model, q=q, point_loads=tuple((int(m), float(a), float(P)) for m, a, P in point_loads),
                   nodal_loads=tuple((int(d), float(v)) for d, v in nodal_loads))
```

- [ ] **Step 4: Write `src/crackfem/gridsolve.py`**

```python
"""Direct stiffness assembly and partition solve of a grid (course method, Logan §5.4).

Member end forces are local: f^ = k^ T d - f^_eq = [V1, T1, M1, V2, T2, M2] acting on the member
(Logan's f^1y, m^1x, m^1z, f^2y, m^2x, m^2z). The internal sagging moment follows the beam rule
of solve.moment_at, M(x) = -M1 + V1 x + q x^2/2 + sum P (x - a), and the internal torque is -T1
(Logan Eq. 5.4.9)."""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from .elements import grid_k_local, grid_point_local, grid_T, grid_udl_local
from .gridmodel import GridModel


@dataclass(frozen=True)
class GridStiffness:
    K: np.ndarray            # assembled global stiffness
    kT: np.ndarray           # (n_members, 6, 6): k^ T of each member (local end forces = kT d - f^_eq)
    free: np.ndarray
    fixed: np.ndarray
    cho: tuple               # Cholesky factor of K_ff (K_ff is symmetric positive definite)


@dataclass(frozen=True)
class GridSolution:
    D: np.ndarray            # nodal displacements (w, theta_x, theta_z per node)
    R: np.ndarray            # reactions at restrained DOFs, zero elsewhere
    end_forces: np.ndarray   # (n_members, 6): [V1, T1, M1, V2, T2, M2], local


def member_dofs(model: GridModel) -> np.ndarray:
    i, j = model.members[:, 0], model.members[:, 1]
    return np.stack([3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2], axis=1)


def _per_member(values, model: GridModel, name: str) -> np.ndarray:
    v = np.broadcast_to(np.asarray(values, dtype=float), (model.n_members,))
    if not np.all(np.isfinite(v)) or np.any(v <= 0.0):
        raise ValueError(f"{name} must be finite and positive")
    return v


def stiffness(model: GridModel, EI, GJ) -> GridStiffness:
    EI, GJ = _per_member(EI, model, "EI"), _per_member(GJ, model, "GJ")
    L = model.lengths
    C, S = model.cosines
    dofs = member_dofs(model)
    K = np.zeros((model.n_dof, model.n_dof))
    kT = np.empty((model.n_members, 6, 6))
    for m in range(model.n_members):
        T = grid_T(C[m], S[m])
        kT[m] = grid_k_local(EI[m], GJ[m], L[m]) @ T
        K[np.ix_(dofs[m], dofs[m])] += T.T @ kT[m]            # Logan Eq. 5.4.19: k = T^T k^ T
    fixed = np.array(model.fixed_dofs, dtype=int)
    free = np.setdiff1d(np.arange(model.n_dof), fixed)
    return GridStiffness(K, kT, free, fixed, cho_factor(K[np.ix_(free, free)]))


def local_loads(model: GridModel) -> np.ndarray:
    """(n_members, 6) consistent local loads of each member's UDL and point loads."""
    L = model.lengths
    q = model.udl()
    f = np.array([grid_udl_local(q[m], L[m]) for m in range(model.n_members)])
    for m, a, P in model.point_loads:
        f[m] += grid_point_local(P, a, L[m])
    return f


def load_vector(model: GridModel) -> tuple[np.ndarray, np.ndarray]:
    """Global load vector F (consistent member loads rotated to global, plus nodal loads) and
    the local member loads."""
    C, S = model.cosines
    fl = local_loads(model)
    dofs = member_dofs(model)
    F = np.zeros(model.n_dof)
    for m in range(model.n_members):
        F[dofs[m]] += grid_T(C[m], S[m]).T @ fl[m]
    for dof, v in model.nodal_loads:
        F[dof] += v
    return F, fl


def solve(model: GridModel, EI=None, GJ=None, stiff: GridStiffness | None = None) -> GridSolution:
    """Partition solve K_ff D_f = F_f; pass a GridStiffness to reuse its factorisation."""
    st = stiffness(model, EI, GJ) if stiff is None else stiff
    F, fl = load_vector(model)
    D = np.zeros(model.n_dof)
    D[st.free] = cho_solve(st.cho, F[st.free])
    R = np.zeros(model.n_dof)
    R[st.fixed] = st.K[st.fixed] @ D - F[st.fixed]
    ef = np.einsum("mij,mj->mi", st.kT, D[member_dofs(model)]) - fl
    return GridSolution(D, R, ef)


def moment_at(model: GridModel, sol: GridSolution, m: int, xl: float) -> float:
    """Internal sagging moment at distance xl from node i of member m."""
    V1, M1 = sol.end_forces[m, 0], sol.end_forces[m, 2]
    out = -M1 + V1 * xl + model.udl()[m] * xl * xl / 2.0
    for mm, a, P in model.point_loads:
        if mm == m and a < xl:
            out += P * (xl - a)
    return float(out)


def torque(sol: GridSolution, m: int) -> float:
    """Internal torque of member m (constant along it), Logan Eq. 5.4.9."""
    return float(-sol.end_forces[m, 1])


def midpoint_moments(model: GridModel, sol: GridSolution, members) -> np.ndarray:
    L = model.lengths
    return np.array([moment_at(model, sol, int(m), L[int(m)] / 2.0) for m in members])


def line_node_moments(model: GridModel, sol: GridSolution, members) -> np.ndarray:
    """Sagging moment at every node of a chain of members (member k runs from node k to k + 1).
    At an interior node the moment jumps by the moment the cross members deliver there, so the
    node value is the mean of the two sides; the end nodes take their one side."""
    members = [int(m) for m in members]
    L = model.lengths
    left = [moment_at(model, sol, m, 0.0) for m in members]
    right = [moment_at(model, sol, m, L[m]) for m in members]
    inner = [(r + l) / 2.0 for r, l in zip(right[:-1], left[1:])]
    return np.array([left[0], *inner, right[-1]])
```

- [ ] **Step 5: Run the tests**

Run: `python -m pytest tests/test_gridsolve.py tests/test_grid_element.py -q`
Expected: `18 passed` (7 + 11).

- [ ] **Step 6: Commit**

```bash
git add src/crackfem/gridmodel.py src/crackfem/gridsolve.py tests/test_gridsolve.py
git commit -m "feat: grid assembly and partition solve; Logan Example 5.5 through the full path" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 4: Generic secant loop that reports non-convergence

**Files:**
- Modify: `src/crackfem/nonlinear.py` (append; `solve_secant` of Stage 1 is unchanged)
- Test: `tests/test_secant_iterate.py`

`solve_secant` (Stage 1) is tied to `BeamModel` and raises on non-convergence. The grillage needs the same Ch11 direct iteration for any model, with optional under-relaxation, and a result that says whether it converged (spec §6: report non-convergence). The test pins it to `solve_secant` on the Stage-1 beams.

- [ ] **Step 1: Write the failing test** — `tests/test_secant_iterate.py`

```python
"""secant_iterate: the generic secant loop (any FE model) that returns non-convergence instead
of raising. Checked against the Stage-1 solve_secant on the same beams."""
import numpy as np
import pytest

from crackfem.model import BeamModel, simply_supported, with_loads
from crackfem.nonlinear import secant_iterate, solve_secant
from crackfem.solve import midpoint_moments, solve

L, EI0, W = 24.5, 7.5e9, 28e3
MREF = W * L**2 / 8


def softening(M):
    return EI0 / (1.0 + np.abs(M) / MREF)


def beam_analyse(model):
    def run(EI):
        sol = solve(model, EI)
        return sol, midpoint_moments(model, sol)
    return run


def test_determinate_beam_one_pass():
    m = with_loads(simply_supported(L, 20), q=-W)
    res = secant_iterate(beam_analyse(m), softening, np.full(m.n_el, EI0))
    assert res.converged and res.reason == "converged"
    assert res.iterations == 1


def test_matches_solve_secant_on_an_indeterminate_beam():
    n = 20
    m = with_loads(BeamModel(x=np.linspace(0, L, n + 1), fixed_dofs=(0, 1, 2 * n)), q=-W)
    ref = solve_secant(m, softening, EI0, tol=1e-10)
    res = secant_iterate(beam_analyse(m), softening, np.full(n, EI0), tol=1e-10)
    assert res.converged and res.iterations == ref.iterations >= 2
    np.testing.assert_allclose(res.EI, ref.EI, rtol=1e-12)
    np.testing.assert_allclose(res.history, ref.history, rtol=1e-9, atol=1e-15)
    assert res.ei_residual[-1] <= 1e-6


def test_two_cycle_is_returned_as_not_converged_with_its_history():
    n = 40
    m = with_loads(BeamModel(x=np.linspace(0, 2 * L, n + 1), fixed_dofs=(0, n, 2 * n)), q=-W)
    Mc = 0.75 * np.max(np.abs(midpoint_moments(m, solve(m, EI0))))
    step = lambda M: np.where(np.abs(M) < Mc, EI0, EI0 / 4)
    res = secant_iterate(beam_analyse(m), step, np.full(n, EI0), max_iter=30)
    assert not res.converged and res.reason == "max_iter"
    assert res.iterations == 30 and len(res.history) == 30 == len(res.ei_residual)


def test_under_relaxation_reaches_the_same_fixed_point():
    n = 20
    m = with_loads(BeamModel(x=np.linspace(0, L, n + 1), fixed_dofs=(0, 1, 2 * n)), q=-W)
    a = secant_iterate(beam_analyse(m), softening, np.full(n, EI0), tol=1e-10)
    b = secant_iterate(beam_analyse(m), softening, np.full(n, EI0), tol=1e-10, relax=0.5)
    assert b.converged and b.iterations > a.iterations
    np.testing.assert_allclose(b.EI, a.EI, rtol=1e-5)          # both stop at tol_EI = 1e-6
    with pytest.raises(ValueError):
        secant_iterate(beam_analyse(m), softening, np.full(n, EI0), relax=0.0)
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_secant_iterate.py -q`
Expected: FAIL — `ImportError: cannot import name 'secant_iterate' from 'crackfem.nonlinear'`.

- [ ] **Step 3: Append to `src/crackfem/nonlinear.py`**

Add at the end of the file (the module already imports `dataclass`, `Callable` and `numpy as np`):

```python
@dataclass(frozen=True)
class SecantResult:
    converged: bool
    solution: object        # whatever analyse returned for the final stiffness
    EI: np.ndarray          # stiffness used in the final solve
    M: np.ndarray           # moments of the final solve (what stiffness_fn reads)
    iterations: int         # number of stiffness updates
    history: tuple          # relative moment change after each update
    ei_residual: tuple      # max|stiffness_fn(M) - EI| / max(EI) after each update
    reason: str             # "converged" | "max_iter"


def secant_iterate(analyse: Callable[[np.ndarray], tuple], stiffness_fn: Callable[[np.ndarray], np.ndarray],
                   EI_start, tol: float = 1e-6, tol_EI: float = 1e-6, max_iter: int = 100,
                   relax: float = 1.0) -> SecantResult:
    """Direct iteration with secant stiffness for any FE model (Ch11): analyse(EI) returns
    (solution, M) and stiffness_fn(M) the secant EI belonging to M. Each update uses
    EI <- relax stiffness_fn(M) + (1 - relax) EI (relax = 1: plain secant; < 1: under-relaxed).
    Stops when the moments stop changing (relative change <= tol) and the stiffness belonging to
    the final moments equals the stiffness used (residual <= tol_EI). Never raises for
    non-convergence: returns converged=False with the whole history, so the caller can report it."""
    if not 0.0 < relax <= 1.0:
        raise ValueError("relax must be in (0, 1]")
    EI = np.asarray(EI_start, dtype=float).copy()
    sol, M = analyse(EI)
    history, resid = [], []
    for it in range(1, max_iter + 1):
        EI = relax * np.asarray(stiffness_fn(M), dtype=float) + (1.0 - relax) * EI
        sol, M_new = analyse(EI)
        change = float(np.max(np.abs(M_new - M)) / max(float(np.max(np.abs(M_new))), 1e-12))
        r = float(np.max(np.abs(np.asarray(stiffness_fn(M_new), dtype=float) - EI)) / np.max(EI))
        history.append(change)
        resid.append(r)
        M = M_new
        if change <= tol and r <= tol_EI:
            return SecantResult(True, sol, EI, M, it, tuple(history), tuple(resid), "converged")
    return SecantResult(False, sol, EI, M, max_iter, tuple(history), tuple(resid), "max_iter")
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_secant_iterate.py tests/test_nonlinear.py tests/test_stage1_golden.py -q`
Expected: `14 passed` (4 + 5 + 5).

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/nonlinear.py tests/test_secant_iterate.py
git commit -m "feat: generic secant iteration with under-relaxation and a non-converged result" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

**Group B4-a ends here.** Review: element signs against the Conventions; the Logan numbers come from the fixture (printed values), not from the code; `python -m pytest -q` → all passed (baseline + 33).

---

## Task 5: S01 deck as a grillage (layout, Hambly properties, loads, truck placements)

**Files:**
- Create: `src/crackfem/deck.py`
- Test: `tests/test_deck.py`, `tests/test_grillage_linear.py`

`test_deck.py` checks the layout, the closed-form properties and the statics of every load mapping (force and moment about x kept exactly). `test_grillage_linear.py` holds the spec's independent linear checks: identical girders under equal load each carry qL²/8 and deflect 5qL⁴/384EI (K symmetric, ΣR = total load); a lone girder line reproduces the Stage-1 `Girder1D` (linear and through its own secant stiffness, one pass); a torsion-free deck reproduces Courbon's rule with rigid cross members and the lever rule with flexible ones, and the real deck lies between them.

Probe values (n_el = 10): wearing surface + barriers give 6.64 kN/m on each exterior girder (plus a torque of 3.83 kN·m/m) and 2.42 kN/m on each interior girder; roadway 10.02 m → 2 design lanes; 6757 transverse placements at a 0.05 m step.

- [ ] **Step 1: Write the failing tests** — `tests/test_deck.py`

```python
"""Deck layout, Hambly properties and load mapping: counts, closed-form values and statics."""
import numpy as np
import pytest

from crackfem.deck import (area_to_girders, fixed_properties, girder_torsion_constant, layout,
                           lever_weights, line_loads_model, line_to_girders, n_lanes,
                           stage_b_dead_lines, torsion_constant_open, truck_placements, wheel_loads)
from crackfem.gridsolve import load_vector
from crackfem.materials import ec2_ecm
from crackfem.params import FT, INCH, Inputs

INP = Inputs()


@pytest.fixture(scope="module")
def deck():
    return layout(INP, 10)


def test_layout_counts_and_supports(deck):
    ng, n = INP.n_girders, 10
    assert deck.model.n_nodes == ng * (n + 1)
    assert deck.model.n_members == ng * n + (n + 1) * (ng - 1)
    assert deck.kind.count("diaphragm") == 3 * (ng - 1)
    assert len(deck.model.fixed_dofs) == 2 * ng and all(d % 3 == 0 for d in deck.model.fixed_dofs)
    assert deck.strip_width[deck.cross_members].sum() == pytest.approx(INP.span * (ng - 1))
    np.testing.assert_allclose(deck.girder_widths(), [1.1 + 3 * FT, 2.2, 2.2, 2.2, 1.1 + 3 * FT])
    with pytest.raises(ValueError):
        layout(INP, 7)


def test_type_iv_torsion_constant_is_sum_b_t3_over_3():
    assert girder_torsion_constant(INP) == pytest.approx((26 + 20 + 38) * 8**3 / 3 * INCH**4, rel=1e-12)
    assert torsion_constant_open([(0.2, 1.0)]) == pytest.approx(1.0 * 0.2**3 / 3)


def test_hambly_slab_strip_and_girder_torsion(deck):
    EI, GJ = fixed_properties(deck)
    Es = ec2_ecm(INP.fck_slab)
    Gs = Es / 2.4
    t, b = INP.slab_thickness, INP.span / 10
    m = int(deck.cross_members[3, 1])                       # an interior slab strip
    assert deck.kind[m] == "slab"
    assert EI[m] == pytest.approx(Es * b * t**3 / 12) and GJ[m] == pytest.approx(Gs * b * t**3 / 6)
    d = int(deck.cross_members[5, 1])                       # the midspan diaphragm
    assert deck.kind[d] == "diaphragm"
    assert EI[d] == pytest.approx(Es * (b * t**3 + 0.25 * 1.10**3) / 12)
    assert np.all(np.isnan(EI[deck.girder_members.ravel()]))
    Gg = ec2_ecm(INP.fck_girder) / 2.4
    assert GJ[deck.girder_members[2, 0]] == pytest.approx(Gg * girder_torsion_constant(INP) + Gs * 2.2 * t**3 / 6)


def test_lever_rule_and_overhang_torque():
    z = np.array([0.0, 2.2, 4.4])
    assert lever_weights(z, 0.55) == [(0, 0.75, 0.0), (1, 0.25, 0.0)]
    assert lever_weights(z, -0.5) == [(0, 1.0, -0.5)]
    q, t = line_to_girders(z, -0.5, 10.0)
    assert q.tolist() == [10.0, 0.0, 0.0] and t.tolist() == [-5.0, 0.0, 0.0]


@pytest.mark.parametrize("za, zb", [(-0.6, 9.4), (0.3, 1.9), (-0.9, 0.0), (1.0, 7.0)])
def test_area_load_mapping_keeps_force_and_moment(za, zb):
    z = INP.girder_spacing * np.arange(5)
    w = 1.1e3
    q, t = area_to_girders(z, za, zb, w)
    assert q.sum() == pytest.approx(w * (zb - za), rel=1e-12)
    # moment about the x axis through z = 0: m_x = z P for a downward load P
    assert (q * z).sum() + t.sum() == pytest.approx(w * (zb**2 - za**2) / 2, rel=1e-12, abs=1e-9)


def test_stage_b_dead_load_totals(deck):
    q, t = stage_b_dead_lines(deck)
    lo, hi = deck.curbs
    total = INP.unit_weight_wearing * INP.wearing_thickness * (hi - lo) + INP.n_barriers * INP.barrier_weight
    assert q.sum() == pytest.approx(total, rel=1e-12)
    np.testing.assert_allclose(q, q[::-1], rtol=1e-12)                    # symmetric deck
    assert t[0] < 0.0 and t[-1] == pytest.approx(-t[0]) and np.all(t[1:-1] == 0.0)
    F, _ = load_vector(line_loads_model(deck, q, t))
    assert -F[0::3].sum() == pytest.approx(total * INP.span, rel=1e-12)
    assert F[1::3].sum() == pytest.approx(t.sum() * INP.span, abs=1e-6)


def test_wheel_loads_keep_force_and_moment(deck):
    axles = [(3.0, 100e3), (12.25, 50e3)]
    lines = (-0.3, 1.2, 5.0)                                              # one wheel line on the overhang
    pts, nodal = wheel_loads(deck, axles, lines, 1.5)
    total = 1.5 * sum(P for _, P in axles) / 2 * len(lines)
    assert -sum(P for _, _, P in pts) == pytest.approx(total, rel=1e-12)
    z_of = {int(m): deck.z[g] for g in range(deck.n_girders) for m in deck.girder_members[g]}
    moment = sum(-P * z_of[m] for m, _, P in pts) + sum(v for dof, v in nodal if dof % 3 == 1)
    assert moment == pytest.approx(1.5 * sum(P for _, P in axles) / 2 * sum(lines), rel=1e-12)


def test_truck_placements_respect_curbs_lanes_and_clearance(deck):
    lo, hi = deck.curbs
    assert n_lanes(deck) == 2                                              # 10.0 m roadway
    ps = truck_placements(deck, step=0.1)
    edge, sw, cw = 2 * FT, 6 * FT, 10 * FT
    assert {len(p) for p in ps} == {2, 4}
    for p in ps:
        assert min(p) >= lo + edge - 1e-9 and max(p) <= hi - edge + 1e-9
        for a, b in zip(p[0::2], p[1::2]):
            assert b - a == pytest.approx(sw)
        if len(p) == 4:
            assert p[2] - p[1] >= 2 * edge - 1e-9                          # clearance widths do not overlap
    has = lambda target: any(len(p) == len(target) and np.allclose(p, target, atol=1e-9) for p in ps)
    assert has((lo + edge, lo + edge + sw)) and has((hi - edge - sw, hi - edge))   # both curb positions
    mirrored = {tuple(round(lo + hi - v, 6) for v in reversed(p)) for p in ps}
    assert mirrored == {tuple(round(v, 6) for v in p) for p in ps}
    assert cw == pytest.approx(2 * edge + sw)
```

and `tests/test_grillage_linear.py`:

```python
"""Linear grillage checks against independent closed forms (spec §6, grillage rows):
identical girders share equally, a lone girder line is the Stage-1 girder, and a torsion-free
deck between flexible and rigid cross members reproduces the lever rule and Courbon's rule."""
import numpy as np
import pytest

from crackfem.deck import fixed_properties, layout
from crackfem.girder1d import Girder1D
from crackfem.gridmodel import GridModel, with_loads
from crackfem.gridsolve import line_node_moments, midpoint_moments, solve, stiffness
from crackfem.nonlinear import secant_iterate
from crackfem.params import Inputs

INP = Inputs()
N = 10
EI_G = 2.6e10          # any girder EI; the checks are ratios or closed forms in EI


@pytest.fixture(scope="module")
def deck():
    return layout(INP, N)


def _props(deck, cross_scale=1.0, gj_scale=1.0):
    EI, GJ = fixed_properties(deck)
    EI = EI * cross_scale
    EI[deck.girder_members.ravel()] = EI_G
    return EI, GJ * gj_scale


def test_identical_girders_under_equal_load_share_equally(deck):
    q = 20e3
    gm = deck.girder_members.ravel()
    qm = np.zeros(deck.model.n_members)
    qm[gm] = -q
    model = with_loads(deck.model, q=qm)
    EI, GJ = _props(deck)
    st = stiffness(model, EI, GJ)
    np.testing.assert_allclose(st.K, st.K.T, rtol=1e-10, atol=1e-2)
    sol = solve(model, stiff=st)
    L = INP.span
    for g in range(deck.n_girders):
        m = line_node_moments(model, sol, deck.girder_members[g])[deck.mid]
        assert m == pytest.approx(q * L**2 / 8, rel=1e-9)
        assert -sol.D[3 * deck.node(g, deck.mid)] == pytest.approx(5 * q * L**4 / (384 * EI_G), rel=1e-9)
    assert sol.R[0::3].sum() == pytest.approx(deck.n_girders * q * L, rel=1e-12)


def test_a_lone_girder_line_is_the_stage1_girder():
    g = Girder1D(n_el=N)
    x = g.model_service.x
    line = GridModel(np.stack([x, np.zeros_like(x)], axis=1), np.array([(k, k + 1) for k in range(N)]),
                     (0, 1, 3 * N))
    dx = x[1] - x[0]
    pts = [(min(int(xp // dx), N - 1), xp - dx * min(int(xp // dx), N - 1), P) for xp, P in g.model_service.point_loads]
    model = with_loads(line, q=-g.q_b, point_loads=pts)
    GJ = 1.0e9
    for eta in (1.0, 0.36):
        ref = g.analyse(eta, g.model_service)

        def run(EI):
            sol = solve(model, EI, GJ)
            return sol, midpoint_moments(model, sol, range(N))

        res = secant_iterate(run, g.stiffness_fn(eta, g.precracked_zone(eta)), np.full(N, g.EI_I))
        assert res.converged and res.iterations == ref.iterations == 1
        np.testing.assert_allclose(res.EI, ref.EI, rtol=1e-9)
        np.testing.assert_allclose(res.solution.D[0::3], ref.solution.D[0::2], rtol=1e-9, atol=1e-12)


def _girder_shares(deck, model, EI, GJ):
    sol = solve(model, EI, GJ)
    m = np.array([line_node_moments(model, sol, deck.girder_members[g])[deck.mid] for g in range(deck.n_girders)])
    return m / m.sum()


def test_torsion_free_deck_between_lever_rule_and_courbon(deck):
    P = 100e3
    model = with_loads(deck.model, nodal_loads=[(3 * deck.node(0, deck.mid), -P)])   # on G1 at midspan
    zbar = deck.z - deck.z.mean()
    courbon = 1 / deck.n_girders + zbar[0] * zbar / (zbar**2).sum()
    rigid = _girder_shares(deck, model, *_props(deck, cross_scale=1e4, gj_scale=1e-6))
    np.testing.assert_allclose(rigid, courbon, atol=5e-3)
    flexible = _girder_shares(deck, model, *_props(deck, cross_scale=1e-8, gj_scale=1e-6))
    assert flexible[0] == pytest.approx(1.0, abs=1e-3)                      # lever rule: all on G1
    real = _girder_shares(deck, model, *_props(deck))
    assert courbon[0] < real[0] < flexible[0]                               # a real deck lies between
```

- [ ] **Step 2: Run them to verify they fail**

Run: `python -m pytest tests/test_deck.py tests/test_grillage_linear.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'crackfem.deck'`.

- [ ] **Step 3: Write `src/crackfem/deck.py`**

```python
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


def fixed_properties(deck: Deck) -> tuple[np.ndarray, np.ndarray]:
    """(EI, GJ) per member. Cross members are linear: slab strip I = b t^3/12, J = b t^3/6
    (Hambly); a diaphragm adds its web (T section of slab concrete, J adds h b^3/3). Girder members
    get GJ = G_girder J_TypeIV + G_slab b t^3/6 with b = their slab width; their EI is NaN here
    (it comes from the section models). A diaphragm is its web (depth x thickness, slab
    concrete) added to the slab strip of its station: not composite with the slab, so its
    properties do not depend on the mesh (a lower bound on its stiffness)."""
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
    for m in deck.cross_members.ravel():
        b = deck.strip_width[m]
        EI[m] = Es * b * t**3 / 12.0
        GJ[m] = Gs * b * t**3 / 6.0
        if deck.kind[m] == "diaphragm":         # the web alone, added to the strip of its station
            EI[m] += Es * bd * hd**3 / 12.0
            GJ[m] += Gs * torsion_constant_open([(hd, bd)])
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


def truck_placements(deck: Deck, step: float = 0.05) -> list[tuple[float, ...]]:
    """Transverse truck placements as tuples of wheel-line z: one truck, and (if the roadway has
    two or more design lanes) two trucks side by side. A truck occupies its clearance width
    (10 ft) with its wheel lines 2 ft inside it, the clearance widths stay between the curbs and
    do not overlap (AASHTO Standard truck in its lane). Edge positions are always included."""
    inp = deck.inp
    lo, hi = deck.curbs
    cw = inp.truck_clearance_ft * FT
    edge = inp.wheel_curb_clearance_ft * FT
    sw = inp.wheel_spacing_ft * FT
    if hi - lo < cw:
        raise ValueError("the roadway is narrower than one truck")
    starts = np.unique(np.concatenate([np.arange(lo, hi - cw, step), [hi - cw]]))
    lines = lambda s: (s + edge, s + edge + sw)
    out = [lines(s) for s in starts]
    if n_lanes(deck) >= 2 and hi - lo >= 2 * cw:
        for s1 in starts:
            s2s = starts[starts >= s1 + cw - 1e-12]
            out += [lines(s1) + lines(s2) for s2 in s2s]
    # add the mirror image of every placement, so that a symmetric deck gets symmetric answers
    mirror = [tuple(sorted(lo + hi - z for z in p)) for p in out]
    rounded = {tuple(round(v, 9) for v in p): p for p in out + mirror}
    return [tuple(float(v) for v in rounded[key]) for key in sorted(rounded)]
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_deck.py tests/test_grillage_linear.py -q`
Expected: `14 passed` (11 + 3).

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/deck.py tests/test_deck.py tests/test_grillage_linear.py
git commit -m "feat: S01 deck grillage (Hambly properties, lever-rule loads, truck placements) with linear checks" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 6: `Grillage` — per-girder Stage-1 sections, strand design, critical placements, envelope, secant analysis

**Files:**
- Create: `src/crackfem/grillage.py`
- Test: `tests/test_grillage.py`

`Girder1D(replace(inp, girder_spacing=slab_width), n_el, n_strands=count)` is the Stage-1 girder of each line (exterior slab width S/2 + overhang = 2.01 m, so its stage-A load and composite section are the exterior girder's own). The strand count of every girder comes from the Stage-1 zero-tension rule (`design_prestress_staged`) under that girder's own grillage service moment at η = 1, iterated until the counts stop changing (the grillage stiffness depends a little on the counts). The truck's transverse placement is chosen per girder (largest midspan stage-B truck moment on the linear intact grillage, by superposition of one solve per girder line) and held fixed for every η.

Probe values (n_el = 20): strands (31, 29, 27, 29, 31) after two rounds; soffit stress under service −0.52 / −0.28 / −0.52 MPa (G1 / G2 / G3); truck midspan moment at the critical placement 1839 / 1433 / 1088 kN·m; G1's critical placement is two trucks against the left curb.

- [ ] **Step 1: Write the failing test** — `tests/test_grillage.py`

```python
"""Grillage of S01 at eta = 1: strand design, truck placement, sanity bands, linearity."""
import numpy as np
import pytest

from crackfem.deck import lever_weights
from crackfem.grillage import Grillage
from crackfem.gridsolve import line_node_moments, midpoint_moments, solve
from crackfem.loads import truck_fraction, truck_position_max_midspan
from crackfem.params import FT
from crackfem.staging import prestress_staged


@pytest.fixture(scope="module")
def G():
    return Grillage(n_el=10)


def test_every_girder_has_no_soffit_tension_at_service_and_the_count_is_minimal(G):
    mid = G.deck.mid
    assert all(s <= 0.0 for s in G.soffit_service)
    for g, gd in enumerate(G.girders):
        M_A = float(gd.M_A_nd[mid])
        fewer = prestress_staged(gd.make_section, gd.is_girder, gd.h_girder, G.n_strands[g] - 1, G.inp.strand_area,
                                 G.fpe, M_A, M_A + G.M_B_service_mid_I[g])
        assert fewer.soffit_stress_service > 0.0
    assert G.design_rounds[-1] == G.n_strands


def test_symmetric_deck_gives_symmetric_design_and_moments(G):
    assert G.n_strands == G.n_strands[::-1]
    np.testing.assert_allclose(G.M_B_service_mid_I, G.M_B_service_mid_I[::-1], rtol=1e-9)
    np.testing.assert_allclose(G.envelope_el, G.envelope_el[::-1], rtol=1e-9)


def test_superposed_placement_moments_equal_a_direct_solve(G):
    for p in (0, len(G.placements) // 3, len(G.placements) - 1):
        model = G.truck_model(G.placements[p])
        sol = solve(model, stiff=G.stiff_I)
        direct = [line_node_moments(model, sol, G.deck.girder_members[g])[G.deck.mid] for g in range(G.deck.n_girders)]
        np.testing.assert_allclose(G.M_truck_mid[p], direct, rtol=1e-9)


def test_exterior_girder_is_loaded_from_its_curb(G):
    lo, _ = G.deck.curbs
    lines = G.placements[G.critical[0]]
    assert min(lines) == pytest.approx(lo + 2 * FT)
    assert len(lines) == 4                                                  # two trucks govern G1


def test_truck_share_sanity_bands(G):
    one_truck = truck_position_max_midspan(G.inp)[1] * G.truck_scale
    frac = G.M_truck_mid[[G.critical[g] for g in range(5)], range(5)] / one_truck
    # interior girders: below the AASHTO S/5.5 fraction (a code bound), and not far below it
    assert 0.5 < max(frac[1:4]) / truck_fraction(G.inp) < 1.1
    # exterior girder: between the lever rule (hinged slab) and Courbon (rigid deck), same wheels
    z, lines = G.deck.z, G.placements[G.critical[0]]
    lever = sum(s / 2 for zw in lines for g, s, _ in lever_weights(z, zw) if g == 0)
    zbar = z - z.mean()
    courbon = sum((1 / z.size + (zw - z.mean()) * zbar[0] / (zbar**2).sum()) / 2 for zw in lines)
    assert lever < frac[0] < courbon


def test_intact_analysis_is_the_linear_solve(G):
    model = G.service_model(0)
    an = G.analyse(model, (1.0,) * 5)
    lin = solve(model, G.EI_I, G.GJ0)
    assert an.result.converged and an.result.iterations == 1
    np.testing.assert_allclose(an.result.solution.D, lin.D, rtol=1e-9, atol=1e-12)
    assert not an.precracked.any()


def test_envelope_covers_the_service_moments_and_nothing_is_precracked_at_eta_1(G):
    for g in range(5):
        model = G.truck_model(G.placements[G.critical[g]], base=G.model_dead)
        sol = solve(model, stiff=G.stiff_I)
        service = G.girders[g].M_A_el + midpoint_moments(model, sol, G.deck.girder_members[g])
        assert np.all(G.envelope_el[g] >= service - 1e-6 * np.abs(service).max())
        assert not G.precracked(g, 1.0).any()
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_grillage.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'crackfem.grillage'`.

- [ ] **Step 3: Write `src/crackfem/grillage.py`**

```python
"""Stage 2: the chain on the grillage of the S01 deck (spec §3, Stage 2).

Each girder g is a Stage-1 Girder1D built with its own slab width (exterior: S/2 + overhang;
interior: S) and its own strand count, so its staged sections, stage-A moments (statics on the
bare girder), Mcr(eta) and crack width are exactly the Stage-1 physics (reused, not copied).
Stage B (wearing surface, barriers, truck) acts on the grillage. The secant loop updates the EI
of every girder element from that element's own total moment; cross members and torsion stay
linear. The cracked girder Gk has eta < 1, the others eta = 1 under the same rule.

Pre-cracked elements: an element whose total service moment under the moving-truck envelope
(linear intact grillage, the truck at girder g's critical transverse placement) exceeds Mcr(eta)
starts pre-cracked (fct = 0, as Stage 1); inside an analysis an element that exceeds the sound
Mcr in any iteration stays pre-cracked for the rest of that analysis (sticky flag)."""
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
    result: SecantResult
    etas: tuple
    precracked: np.ndarray        # (n_girders, n_el) after the analysis (sticky flags included)
    M_B_mid: np.ndarray           # stage-B sagging moment at each girder's midspan node
    deflection_mid: np.ndarray    # downward midspan deflection of each girder (m)


class Grillage:
    def __init__(self, inp: Inputs = Inputs(), n_el: int = 20, dy_max: float = 1e-3,
                 n_strands=None, placement_step: float = 0.05, max_design_rounds: int = 6):
        """n_strands=None designs every girder by the Stage-1 zero-tension rule under its own
        grillage service moment (eta = 1); a tuple analyses those as-built counts."""
        self.inp, self.n_el, self.dy_max = inp, n_el, dy_max
        self.deck = layout(inp, n_el)
        self.widths = self.deck.girder_widths()
        self.EI_fixed, self.GJ0 = fixed_properties(self.deck)
        self.J_girder = girder_torsion_constant(inp)
        self.q_dead, self.t_dead = stage_b_dead_lines(self.deck)
        self.model_dead = line_loads_model(self.deck, self.q_dead, self.t_dead)
        self.axles, _ = truck_position_max_midspan(inp)         # Stage-1 longitudinal placement
        self.truck_scale = (1.0 + impact_factor(inp)) * inp.live_load_factor
        self.placements = truck_placements(self.deck, placement_step)
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

    def stiffness_fn(self, etas, pre0):
        """EI of every girder element from its stage-B midpoint moment (deck order g, k); the
        element's total moment is M_A + M_B and its reference state the end of stage A (Stage-1
        rule). Sticky flag: an element whose total moment exceeds the sound Mcr in any call stays
        pre-cracked (decompression Mcr) from then on."""
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
                    if M > sound[i]:
                        pre[g, i] = True
                    m_cr = dec[i] if pre[g, i] else sound[i]
                    out[g * n + i] = ec2_secant_EI(gd.sec_el[i], e0[g], M, m_cr, beta, M_ref)
            return out

        fn.precracked = pre
        return fn

    def analyse(self, model, etas, gj_factor: float = 1.0, relax: float = 1.0, tol: float = 1e-6,
                tol_EI: float = 1e-6, max_iter: int = 100) -> GrillageAnalysis:
        """Secant iteration on the grillage. gj_factor scales the GJ of girder elements that are
        cracked (EI below uncracked) - a sensitivity option; 1.0 keeps torsion linear."""
        etas = tuple(float(e) for e in etas)
        if len(etas) != self.deck.n_girders:
            raise ValueError("etas needs one value per girder")
        pre0 = np.array([self.precracked(g, e) for g, e in enumerate(etas)])
        fn = self.stiffness_fn(etas, pre0)
        ids, EI_I_g = self.girder_ids, self.EI_I[self.girder_ids]

        def run(EI_g):
            EI = self.EI_fixed.copy()
            EI[ids] = EI_g
            GJ = self.GJ0.copy()
            if gj_factor != 1.0:
                cracked = EI_g < EI_I_g * (1.0 - 1e-9)
                GJ[ids[cracked]] *= gj_factor
            sol = solve(model, EI, GJ)
            return sol, midpoint_moments(model, sol, ids)

        res = secant_iterate(run, fn, EI_I_g, tol=tol, tol_EI=tol_EI, max_iter=max_iter, relax=relax)
        d, mid = self.deck, self.deck.mid
        sol = res.solution
        M_mid = np.array([line_node_moments(model, sol, d.girder_members[g])[mid] for g in range(d.n_girders)])
        defl = np.array([-sol.D[3 * d.node(g, mid)] for g in range(d.n_girders)])
        return GrillageAnalysis(res, etas, fn.precracked.copy(), M_mid, defl)
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_grillage.py tests/test_stage1_golden.py -q`
Expected: `12 passed` (7 + 5), in about 5 s.

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/grillage.py tests/test_grillage.py
git commit -m "feat: grillage of S01 - per-girder Stage-1 sections, strand design, critical truck placement, envelope" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 7: Nonlinear grillage — w_pred(η) and the indicator

**Files:**
- Modify: `src/crackfem/grillage.py` (append methods to class `Grillage`)
- Test: `tests/test_grillage_nonlinear.py`

This is spec P5's point: on the indeterminate grillage the secant loop needs more than one pass, ends with the stiffness that belongs to its own moments, and a cracked girder sheds moment that its neighbours pick up (the sum of the girders' midspan moments stays about the same: load is shed, not lost). Probe (n_el = 20, G1 at η = 0.36): service loop converged in 34 updates; G1's stage-B midspan moment fell from 2115 to 740 kN·m while G2 rose from 1778 to 2868 kN·m; G2 also cracked (3 elements) — the zero-tension design of the neighbours is exceeded when they take G1's share, and the uniform rule handles it.

- [ ] **Step 1: Write the failing test** — `tests/test_grillage_nonlinear.py`

```python
"""The secant loop on the indeterminate grillage (spec P5): more than one pass, consistent EI,
load shed (not lost), sticky pre-cracked flags, w_pred and the inverse through the grillage."""
import numpy as np
import pytest

from crackfem.grillage import Grillage
from crackfem.inverse import estimate_eta, is_monotone_nonincreasing


@pytest.fixture(scope="module")
def G():
    return Grillage(n_el=10)


@pytest.fixture(scope="module")
def damaged(G):
    return G.analyse(G.service_model(0), G.etas(0, 0.3))


def test_loop_needs_more_than_one_pass_and_ends_with_consistent_stiffness(G, damaged):
    r = damaged.result
    assert r.converged and r.iterations >= 2
    assert r.history[-1] <= 1e-6 and r.ei_residual[-1] <= 1e-6
    fn = G.stiffness_fn(G.etas(0, 0.3), damaged.precracked)      # the stiffness belonging to the final moments
    np.testing.assert_allclose(fn(r.M), r.EI, rtol=1e-5)


def test_cracked_girder_sheds_moment_that_the_others_pick_up(G, damaged):
    intact = G.analyse(G.service_model(0), (1.0,) * 5)
    assert damaged.M_B_mid[0] < intact.M_B_mid[0]
    assert damaged.M_B_mid[1] > intact.M_B_mid[1]
    # load is shed, not lost: the girders' midspan moments still add up to about the same
    assert damaged.M_B_mid.sum() == pytest.approx(intact.M_B_mid.sum(), rel=1e-2)
    assert damaged.deflection_mid[0] > intact.deflection_mid[0]


def test_sticky_flag_keeps_an_element_precracked(G):
    gd = G.girders[2]
    etas = G.etas(2, 0.45)
    sound, dec = G._mcr(2, 0.45)
    i = G.deck.mid                                         # element just right of midspan
    n = G.n_el
    assert dec[i] < sound[i]
    M_b = np.zeros(5 * n)
    fn = G.stiffness_fn(etas, np.zeros((5, n), dtype=bool))
    M_b[2 * n + i] = 1.05 * sound[i] - gd.M_A_el[i]        # above the sound Mcr: flag set
    fn(M_b)
    M_b[2 * n + i] = 0.5 * (sound[i] + dec[i]) - gd.M_A_el[i]   # back between decompression and Mcr
    sticky = fn(M_b)[2 * n + i]
    fresh = G.stiffness_fn(etas, np.zeros((5, n), dtype=bool))(M_b)[2 * n + i]
    assert fresh == pytest.approx(gd.EI_I) and sticky < 0.99 * gd.EI_I
    assert fn.precracked[2, i]


def test_w_pred_through_the_grillage(G):
    assert G.w_pred(1.0, 0) == 0.0
    assert G.w_pred(0.2, 0) > 0.0
    assert is_monotone_nonincreasing(G.w_pred_fn(2), np.linspace(0.05, 1.0, 12))
    assert G.w_pred(0.3, 0) == pytest.approx(G.w_pred(0.3, 4), rel=1e-9)    # G1 and G5 are mirror images
    d = G.w_pred_detail(0.3, 0)
    assert d["converged"] and d["iterations"] >= 2


def test_inverse_round_trip_through_the_grillage(G):
    est = estimate_eta(G.w_pred_fn(2), G.w_pred(0.33, 2))
    assert est.status == "ok"
    assert est.eta == pytest.approx(0.33, abs=1e-5)


def test_non_convergence_is_reported_not_raised(G):
    an = G.analyse(G.service_model(0), G.etas(0, 0.3), max_iter=1)
    assert not an.result.converged and an.result.reason == "max_iter" and len(an.result.history) == 1


def test_torsion_sensitivity_and_indicator(G):
    base = G.indicator(0.3, 0)
    soft = G.indicator(0.3, 0, gj_factor=0.5)
    assert base["converged"] == {"dead": True, "service": True} and soft["converged"]["service"]
    assert soft["M_B_service_mid_kNm"]["damaged"][0] != pytest.approx(base["M_B_service_mid_kNm"]["damaged"][0], rel=1e-6)
    assert soft["w_pred_mm"] != pytest.approx(base["w_pred_mm"], rel=1e-6)          # the dead-load solve sees it too
    assert base["cracked_elements"][0] > 0
    assert sum(base["share_service"]["damaged"]) == pytest.approx(1.0)
    assert base["shed_percent_of_Gk"][0] < 0.0 < base["shed_percent_of_Gk"][1]
    assert base["deflection_truck_mm"]["damaged"][0] > base["deflection_truck_mm"]["intact"][0]
    intact = G.indicator(1.0, 0)
    assert intact["cracked_elements"] == [0, 0, 0, 0, 0] and intact["iterations"]["service"] == 1
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_grillage_nonlinear.py -q`
Expected: FAIL — `AttributeError: 'Grillage' object has no attribute 'etas'` (and `w_pred`, `indicator`).

- [ ] **Step 3: Append these methods at the end of class `Grillage` in `src/crackfem/grillage.py`** (indented as class methods, after `analyse`)

```python
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
            "iterations": {"dead": dmg["dead"].result.iterations, "service": s.iterations},
            "history_service": list(s.history), "ei_residual_service": list(s.ei_residual),
            "history_dead": list(dmg["dead"].result.history),
        }
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_grillage_nonlinear.py tests/test_grillage.py -q`
Expected: `14 passed` (7 + 7), in about 15 s.

- [ ] **Step 5: Run the whole suite (B4-b close)**

Run: `python -m pytest -q`
Expected: all passed — baseline + 61 new (`369 passed` with a pipeline configuration, `366 passed, 3 skipped` without).

- [ ] **Step 6: Commit**

```bash
git add src/crackfem/grillage.py tests/test_grillage_nonlinear.py
git commit -m "feat: nonlinear grillage - sticky pre-cracked flags, w_pred through the grillage, load-shed indicator" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

**Group B4-b ends here.** Review: statics of the load mapping (test_deck), the zero-tension design and its minimality (test_grillage), the loop's >1 pass and consistent final stiffness (test_grillage_nonlinear).

---

## Task 8: PyNite cross-check of the linear intact grillage (**gate: approved dependency install**)

**Files:**
- Modify: `requirements.txt`
- Create: `src/crackfem/pynite_check.py`
- Test: `tests/test_pynite_check.py`

PyNiteFEA 2.0.2 is the newest release that does not require numpy ≥ 2.4 (2.1.0 and later declare `numpy>=2.4.0`, which would move the pinned numpy 2.3.5 and with it every Stage-1 number). Its non-optional dependencies are numpy, PrettyTable, scipy and matplotlib; PrettyTable is new to this project. Its `Pynite/__init__.py` imports `pip._vendor.pkg_resources` (present in pip 26.0.1 on the development machine).

- [ ] **Step 1: Install (approved dependency; holds numpy, scipy and matplotlib at their pins)**

Run: `python -m pip install "PyNiteFEA==2.0.2" "numpy==2.3.5" "scipy==1.17.1" "matplotlib==3.11.2"`
then: `python -m pip show PyNiteFEA` and `python -m pip check`
Expected: `Name: PyNiteFEA`, `Version: 2.0.2`; `No broken requirements found.`
If the install fails, or `python -m pip check` reports a conflict, **stop and report BLOCKED** with the full output. Do not try another version or package.

- [ ] **Step 2: Pin it**

Run: `python -m pip freeze | grep -i -E "^(pynitefea|prettytable|wcwidth)=="`
Append every line it prints to `requirements.txt`, as printed (the PyNiteFEA 2.0.2 line and the exact PrettyTable and, if present, wcwidth versions that pip installed; pip may print the names in either case). Then run `python -m pip install -r requirements.txt`
Expected: `Requirement already satisfied` for every line, nothing upgraded or downgraded.

- [ ] **Step 3: Write the failing test** — `tests/test_pynite_check.py`

```python
"""Independent solver: the same linear grid in PyNite (PyNiteFEA, pinned in requirements.txt)."""
from pathlib import Path

import numpy as np
import pytest

from crackfem.gridmodel import GridModel, with_loads
from crackfem.gridsolve import moment_at, solve
from crackfem.grillage import Grillage
from crackfem.pynite_check import pynite_deflections, pynite_grid, pynite_moment

REQ = Path(__file__).resolve().parents[1] / "requirements.txt"


def test_pynite_is_pinned():
    lines = [line.strip().lower() for line in REQ.read_text(encoding="utf-8").splitlines()]
    assert "pynitefea==2.0.2" in lines                         # pip may print the name in either case


def test_simply_supported_member_against_the_closed_form():
    L, EI, q = 8.0, 2.0e7, -10e3
    s = np.linspace(0.0, L, 3)
    model = with_loads(GridModel(np.stack([s, np.zeros(3)], axis=1), np.array([(0, 1), (1, 2)]), (0, 1, 6)), q=q)
    fm = pynite_grid(model, EI, 1.0e7)
    w = pynite_deflections(fm, model)
    assert w[1] == pytest.approx(5 * q * L**4 / (384 * EI), rel=1e-6)
    assert abs(pynite_moment(fm, 1, 0.0)) == pytest.approx(abs(q) * L**2 / 8, rel=1e-6)


@pytest.fixture(scope="module")
def G():
    return Grillage(n_el=10)


def test_intact_grillage_service_case_agrees_with_pynite(G):
    model = G.service_model(0)
    ours = solve(model, G.EI_I, G.GJ0)
    fm = pynite_grid(model, G.EI_I, G.GJ0)
    w = pynite_deflections(fm, model)
    np.testing.assert_allclose(w, ours.D[0::3], rtol=1e-6, atol=1e-6 * np.abs(ours.D[0::3]).max())
    d = G.deck
    right_of_mid = [int(d.girder_members[g, d.mid]) for g in range(d.n_girders)]   # starts at midspan
    ours_m = np.array([moment_at(model, ours, m, 0.0) for m in right_of_mid])
    pyn_m = np.array([pynite_moment(fm, m, 0.0) for m in right_of_mid])
    sign = np.sign(pyn_m[2] / ours_m[2])                         # PyNite's own Mz sign convention
    np.testing.assert_allclose(sign * pyn_m, ours_m, rtol=1e-6, atol=1e-6 * np.abs(ours_m).max())
```

- [ ] **Step 4: Run it to verify it fails**

Run: `python -m pytest tests/test_pynite_check.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'crackfem.pynite_check'`.

- [ ] **Step 5: Write `src/crackfem/pynite_check.py`**

```python
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
```

- [ ] **Step 6: Run the tests**

Run: `python -m pytest tests/test_pynite_check.py -q`
Expected: `3 passed`. The simply supported member is the API smoke test: if it fails, the PyNite call is wrong (read the installed `Pynite/FEModel3D.py` and fix the call, keeping the test). PyNite and this code use the same Euler–Bernoulli element and the same consistent loads, so agreement is expected to round-off; a systematic difference larger than 1e-6 is a finding to investigate and report, not a tolerance to widen.

- [ ] **Step 7: Commit**

```bash
git add requirements.txt src/crackfem/pynite_check.py tests/test_pynite_check.py
git commit -m "feat: PyNite cross-check of the linear grillage (PyNiteFEA 2.0.2 pinned)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

**Group B4-c ends here.** B4 is complete when `python -m pytest -q` passes with 64 new tests (`372 passed` with a pipeline configuration, `369 passed, 3 skipped` without).

---

## Task 9: Tabulated forward map w(η) for the Monte Carlo

**Files:**
- Create: `src/crackfem/wtable.py`
- Test: `tests/test_wtable.py`

Each exact w_pred is a nonlinear grillage solve (0.03–0.3 s in the probe), so 2,000 Monte Carlo inversions per case would need about 25,000 solves. The table evaluates the exact map once on a 99-point grid plus the plateau edge (bisection on the exact map, as Stage 1), interpolates linearly (monotone data stays monotone), and is used as the forward map of the unchanged `uq.monte_carlo_eta`. Probe: largest relative interpolation error at the knot midpoints 0.18 % (G1) and 0.24 % (G3); the Monte Carlo percentiles on the table and the exact inverse at the same widths differ by at most 1.4e-4 in η.

- [ ] **Step 1: Write the failing test** — `tests/test_wtable.py`

```python
"""Tabulated w(eta) for the Monte Carlo: exact at the knots, zero on the plateau, monotone,
and the same Monte Carlo percentiles as the exact map."""
import numpy as np
import pytest

from crackfem.uq import monte_carlo_eta
from crackfem.wtable import build_table, max_rel_error


def w_exact(eta):
    return 1e-3 * max(0.0, 0.6 - eta) ** 1.5


def test_table_is_exact_at_knots_and_zero_on_the_plateau():
    t = build_table(w_exact, np.linspace(0.02, 1.0, 99))
    assert t.eta_dec == pytest.approx(0.6, abs=2e-4)          # widths below 1 nm count as no crack
    for e in t.eta[t.eta < t.eta_dec]:
        assert t(float(e)) == pytest.approx(w_exact(float(e)), rel=1e-12)
    assert t(0.61) == 0.0 and t(1.0) == 0.0
    assert t.is_nonincreasing()


def test_interpolation_error_falls_with_the_grid():
    errs = []
    for n in (12, 25, 99):
        t = build_table(w_exact, np.linspace(0.02, 1.0, n))
        knots = t.eta[t.eta < t.eta_dec]
        errs.append(max_rel_error(t, w_exact, 0.5 * (knots[:-1] + knots[1:])))
    assert errs[0] > errs[1] > errs[2] and errs[2] < 5e-3


def test_monte_carlo_on_the_table_matches_the_exact_map():
    t = build_table(w_exact, np.linspace(0.02, 1.0, 99))
    a = monte_carlo_eta(t, 0.15e-3, 0.04e-3, n=2000).summary()
    b = monte_carlo_eta(w_exact, 0.15e-3, 0.04e-3, n=2000).summary()
    for q in ("p05", "p50", "p95"):
        assert a[q] == pytest.approx(b[q], abs=2e-3)


def test_bad_forward_maps_are_refused():
    with pytest.raises(ValueError, match="not finite"):
        build_table(lambda e: float("nan") if e < 0.3 else w_exact(e), np.linspace(0.02, 1.0, 20))
    with pytest.raises(ValueError, match="ill posed"):
        build_table(lambda e: 1e-3 * max(0.0, 0.6 - e) * (1.0 + 5.0 * (e > 0.3)), np.linspace(0.02, 1.0, 20))
    with pytest.raises(ValueError):
        build_table(w_exact, [0.5, 0.2])
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_wtable.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'crackfem.wtable'`.

- [ ] **Step 3: Write `src/crackfem/wtable.py`**

```python
"""Tabulated forward map w(eta) for the Monte Carlo on the grillage.

Each exact w_pred(eta) is a nonlinear grillage solve, so the Monte Carlo (thousands of inversions)
runs on a table instead: w on a dense eta grid plus the edge of the no-crack plateau (found by
bisection on the exact map), joined by straight lines. Linear interpolation of non-increasing
data is non-increasing, so the inverse stays well posed; above the plateau edge the table returns
exactly zero, as the exact map does. The table is checked against the exact map (max_rel_error)
and the Monte Carlo percentiles are confirmed with the exact inverse (run_stage2)."""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .inverse import eta_decompression


@dataclass(frozen=True)
class WTable:
    eta: np.ndarray        # strictly ascending knots
    w: np.ndarray          # exact w_pred at the knots (m)
    eta_dec: float         # cracked-side edge of the no-crack plateau

    def __call__(self, eta: float) -> float:
        if eta > self.eta_dec:
            return 0.0
        return float(np.interp(eta, self.eta, self.w))

    def is_nonincreasing(self) -> bool:
        return bool(np.all(np.diff(self.w) <= 1e-9 * max(float(self.w.max()), 1e-15)))


def build_table(w_pred: Callable[[float], float], grid, xtol: float = 1e-7) -> WTable:
    """Evaluate w_pred on the grid (ascending, within [grid[0], grid[-1]]) and at the plateau
    edge. Raises ValueError if any width is not finite (e.g. a non-converged secant loop) or the
    widths increase with eta (then the inverse is ill posed and the report must say so)."""
    grid = np.asarray(grid, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.diff(grid) > 0.0):
        raise ValueError("grid must be strictly ascending with at least 2 points")
    edge = eta_decompression(w_pred, float(grid[0]), float(grid[-1]), xtol=xtol)
    knots = np.unique(np.concatenate([grid[grid < edge], [edge]]))
    w = np.array([w_pred(float(e)) for e in knots])
    if not np.all(np.isfinite(w)):
        raise ValueError(f"w_pred is not finite at eta = {knots[~np.isfinite(w)].tolist()}")
    if knots[-1] < grid[-1]:                      # plateau: exact zeros from just above the edge
        knots = np.concatenate([knots, [min(edge + xtol, grid[-1])]])
        w = np.concatenate([w, [0.0]])
    t = WTable(knots, w, float(edge))
    if not t.is_nonincreasing():
        raise ValueError("w_pred increases with eta on the grid: the inverse is ill posed")
    return t


def max_rel_error(table: WTable, w_pred: Callable[[float], float], etas) -> float:
    """Largest |table - exact| / max(exact widths) over test points (e.g. the grid midpoints)."""
    exact = np.array([w_pred(float(e)) for e in etas])
    approx = np.array([table(float(e)) for e in etas])
    return float(np.max(np.abs(approx - exact)) / max(float(np.max(np.abs(exact))), 1e-15))
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_wtable.py tests/test_uq.py tests/test_inverse.py -q`
Expected: all passed (4 new).

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/wtable.py tests/test_wtable.py
git commit -m "feat: tabulated w(eta) with plateau edge for the grillage Monte Carlo" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 10: Stage-2 figures

**Files:**
- Modify: `src/crackfem/figures.py` (append; Stage-1 functions unchanged)
- Test: `tests/test_figures_stage2.py`

- [ ] **Step 1: Write the failing test** — `tests/test_figures_stage2.py`

```python
from crackfem import figures as F


def test_every_stage2_figure_function_writes_a_png(tmp_path):
    band = {"p05": 0.3, "p50": 0.35, "p95": 0.4}
    rows = [{"name": "c1", "eta_true": 0.36, "bands": {"Stage 1": band, "grillage G1": band, "grillage G3": None}},
            {"name": "c2", "eta_true": 0.21, "bands": {"Stage 1": band, "grillage G1": band, "grillage G3": band}}]
    calls = [
        lambda p: F.fig_w_eta_models([0.1, 0.5, 1.0], {"a": [0.4, 0.1, 0.0], "b": [0.3, 0.0, 0.0]}, p,
                                     eta_dec={"a": 0.6}),
        lambda p: F.fig_eta_bands(rows, p, title="t"),
        lambda p: F.fig_girder_bars(["G1", "G2"], {"intact": [1.0, 2.0], "damaged": [0.5, 2.6]}, p, "M", title="t"),
        lambda p: F.fig_history({"G1": [0.5, 0.1, 1e-4, 1e-7], "G3": [0.7, 0.0]}, p, title="t"),
    ]
    for i, draw in enumerate(calls):
        p = tmp_path / f"s{i}.png"
        draw(p)
        assert p.stat().st_size > 1000, i
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_figures_stage2.py -q`
Expected: FAIL — `AttributeError: module 'crackfem.figures' has no attribute 'fig_w_eta_models'`.

- [ ] **Step 3: Append to `src/crackfem/figures.py`**

```python
# --- Stage 2 (grillage) ------------------------------------------------------------------
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
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_eta_bands(rows, path, title=None) -> None:
    """rows: [{"name", "eta_true", "bands": {model label: {"p05", "p50", "p95"} or None}}]: the
    Monte Carlo 5-50-95 % of eta per measured case and model, side by side."""
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    labels = sorted({m for r in rows for m in r["bands"]}, key=lambda s: (s != "Stage 1", s))
    width = 0.8 / max(len(labels), 1)
    handles = {}
    for i, r in enumerate(rows):
        handles["η true (Stage-1 synthetic)"] = ax.plot([i], [r["eta_true"]], "kx", ms=7, mew=1.5, zorder=5)[0]
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


def fig_girder_bars(girders, series, path, ylabel, title=None) -> None:
    """Grouped bars per girder (e.g. intact vs damaged midspan moment or deflection)."""
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    n = len(series)
    width = 0.8 / max(n, 1)
    for j, (label, vals) in enumerate(series.items()):
        xs = [i - 0.4 + width * (j + 0.5) for i in range(len(girders))]
        ax.bar(xs, vals, width=width, color=f"C{j}", label=label)
    ax.axhline(0.0, color="k", lw=0.8)
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
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_figures_stage2.py tests/test_figures.py -q`
Expected: `2 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/figures.py tests/test_figures_stage2.py
git commit -m "feat: Stage-2 figures (w(eta) by model, eta bands, girder bars, secant history)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 11: `scripts/run_stage2.py`

**Files:**
- Create: `scripts/run_stage2.py`
- Test: `tests/test_run_stage2.py`

For each cracked girder (G1, G3) the script tabulates w(η), runs every committed thesis-pipeline case (exact point estimate; Monte Carlo on the table; the table's p05/p50/p95 confirmed by the exact inverse at the same widths; a refusal gives the nominal 3 px gate bound), picks the reference case as Stage 1 does (finest GSD, least damage) and evaluates the indicator at its Monte Carlo median with a deflection band from p95/p05, plus a torsion sensitivity run (GJ of cracked girder elements × 0.5). It adds the PyNite cross-check, the mesh study (n_el = 10, 20, 40, 80 at the designed strand counts) and timing. Stage-1 bands are read from `results/stage1.json` (read only).

- [ ] **Step 1: Write the failing test** — `tests/test_run_stage2.py`

```python
"""run_stage2.run_case: measured, refused and unresolved cases on a tabulated forward map."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from crackfem.wtable import build_table

_spec = importlib.util.spec_from_file_location(
    "run_stage2", Path(__file__).resolve().parents[1] / "scripts" / "run_stage2.py")
run_stage2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_stage2)


def w_exact(eta):
    return 1e-3 * max(0.0, 0.6 - eta)


@pytest.fixture(scope="module")
def table():
    return build_table(w_exact, np.linspace(0.02, 1.0, 99))


def _rec(cracks, gsd_mm=0.02):
    return {"case": {"name": "t", "gsd_mm": gsd_mm, "eta_true": 0.36, "w_true_mm": 0.24, "width_px": 12},
            "parsed": {"gsd_mm": gsd_mm, "min_measurable_mm": 3 * gsd_mm, "cracks": cracks}}


def _crack(gate, w_mm=None, u_mm=None):
    return {"index": 1, "points": 390, "gate": gate, "max_width_px": None if w_mm is None else w_mm / 0.02,
            "max_width_mm": w_mm, "uncertainty_mm": u_mm, "rating": 5 if w_mm else None,
            "rating_label": "x" if w_mm else None}


def test_measured_case_band_contains_the_truth_and_the_exact_check_agrees(table):
    out = run_stage2.run_case(table, w_exact, _rec([_crack("ok", 0.24, 0.04)]))
    assert out["status"] == "measured" and out["eta_status"] == "ok"
    assert out["eta_point"] == pytest.approx(0.36, abs=1e-6)
    assert out["eta_mc"]["p05"] < 0.36 < out["eta_mc"]["p95"]
    assert out["eta_mc_exact_max_abs_diff"] < 1e-4


def test_refusal_gives_the_gate_bound_and_unresolved_gives_none(table):
    ref = run_stage2.run_case(table, w_exact, _rec([_crack("below-measurable-resolution")]))
    assert ref["status"] == "refused" and ref["eta_lower_bound_from_gate"] == pytest.approx(0.6 - 0.06, abs=1e-6)
    unr = run_stage2.run_case(table, w_exact, _rec([_crack("no-valid-points")]))
    assert unr["status"] == "unresolved" and unr["eta_lower_bound_from_gate"] is None
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_run_stage2.py -q`
Expected: FAIL — collection error `FileNotFoundError` for `scripts/run_stage2.py`.

- [ ] **Step 3: Write `scripts/run_stage2.py`**

```python
"""Regenerate every Stage-2 (grillage) result and figure from committed inputs.

Reads data/claimbridge/*_spike.json (thesis-pipeline output) and results/stage1.json (Stage-1
bands, read only) and writes results/stage2.json + results/figures/fig9...fig14. Stage-1 outputs
are not touched. Does not need the thesis pipeline; needs PyNiteFEA (requirements.txt)."""
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
    the table's percentiles confirmed by the exact inverse; refusals give the nominal gate bound."""
    c = rec["case"]
    out = {"name": c["name"], "gsd_mm": c["gsd_mm"], "eta_true_stage1": c["eta_true"]}
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


def main(pynite: bool = True) -> None:
    t_all = time.perf_counter()
    res: dict = {"inputs": asdict(Inputs()), "n_el": N_EL}
    t0 = time.perf_counter()
    g = Grillage(n_el=N_EL)
    t_build = time.perf_counter() - t0
    d = g.deck
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
        ref = min(measured, key=lambda c: (c["gsd_mm"], -c["eta_true_stage1"]))
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
            "indicator": dmg,
            "deflection_service_band_mm": {q: band[v]["deflection_service_mm"]["damaged"][k] for q, v in from_eta.items()},
            "deflection_truck_band_mm": {q: band[v]["deflection_truck_mm"]["damaged"][k] for q, v in from_eta.items()},
            "gj_sensitivity": {"gj_factor": GJ_SENSITIVITY,
                               "M_B_service_mid_kNm": sens["M_B_service_mid_kNm"]["damaged"],
                               "deflection_service_mm": sens["deflection_service_mm"]["damaged"],
                               "w_pred_mm": sens["w_pred_mm"]},
        }
    res["stage1_comparison"] = [
        {"name": c["name"], "eta_true": c["eta_true_stage1"],
         "bands": {"Stage 1": s1_cases[c["name"]].get("eta_mc"),
                   **{f"grillage {n}": next(x for x in res[n]["cases"] if x["name"] == c["name"])["eta_mc"] for n in CRACKED}}}
        for c in res["G1"]["cases"] if c["status"] == "measured"]

    if pynite:
        res["pynite"] = pynite_crosscheck(g)

    mesh = {"N": list(MESH_NS)}
    for name, k in CRACKED.items():
        eta = res[name]["eta_ref"]
        w, dfl = [], []
        for n in MESH_NS:
            gm = g if n == N_EL else Grillage(n_el=n, n_strands=g.n_strands)
            w.append(gm.w_pred(eta, k) * 1e3)
            dfl.append(gm.indicator(eta, k)["deflection_service_mm"]["damaged"][k])
        mesh[name] = {"eta": eta, "w_pred_mm": w, "deflection_service_mm": dfl}
    res["mesh"] = mesh

    t0 = time.perf_counter(); solve(g.service_model(0), g.EI_I, g.GJ0); timing["linear_solve_s"] = time.perf_counter() - t0
    t0 = time.perf_counter(); g.analyse(g.service_model(0), g.etas(0, res["G1"]["eta_ref"])); timing["nonlinear_service_G1_s"] = time.perf_counter() - t0
    timing["total_s"] = time.perf_counter() - t_all
    res["timing_s"] = timing

    out = ROOT / "results"
    figs = out / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    (out / "stage2.json").write_text(json.dumps(jsonable(res), indent=2, ensure_ascii=False, allow_nan=False),
                                     encoding="utf-8")

    F.fig_w_eta_models(GRID, {"Stage 1 (interior girder, S/5.5)": res["stage1_w_curve"]["w_mm"],
                              "grillage, G1 cracked": res["G1"]["w_curve"]["w_mm"],
                              "grillage, G3 cracked": res["G3"]["w_curve"]["w_mm"]},
                       figs / "fig9_grillage_w_eta.png",
                       eta_dec={"Stage 1 (interior girder, S/5.5)": stage1["eta_decompression"],
                                "grillage, G1 cracked": res["G1"]["eta_decompression"],
                                "grillage, G3 cracked": res["G3"]["eta_decompression"]})
    F.fig_eta_bands(res["stage1_comparison"], figs / "fig10_eta_bands.png",
                    title="η per measured case: Stage-1 girder vs grillage (G1, G3 cracked)")
    labels = [f"G{i + 1}" for i in range(d.n_girders)]
    for name, fig_m, fig_d in (("G1", "fig11_moment_share_G1.png", "fig12_deflection_G1.png"),
                               ("G3", "fig11_moment_share_G3.png", "fig12_deflection_G3.png")):
        ind, eta = res[name]["indicator"], res[name]["eta_ref"]
        F.fig_girder_bars(labels, {"intact": ind["M_B_service_mid_kNm"]["intact"],
                                   f"{name} at η = {eta:.3f}": ind["M_B_service_mid_kNm"]["damaged"]},
                          figs / fig_m, "stage-B midspan moment (kN·m)\ndead + truck",
                          title=f"{name} cracked: load shed to the other girders")
        F.fig_girder_bars(labels, {"intact": ind["deflection_service_mm"]["intact"],
                                   f"{name} at η = {eta:.3f}": ind["deflection_service_mm"]["damaged"]},
                          figs / fig_d, "midspan deflection (mm, down)\nstage-B dead + truck",
                          title=f"{name} cracked: limit L/{g.inp.deflection_limit_ratio:g} applies to the truck part")
    F.fig_history({f"{n} service (η = {res[n]['eta_ref']:.3f})": res[n]["indicator"]["history_service"] for n in CRACKED}
                  | {f"{n} dead": res[n]["indicator"]["history_dead"] for n in CRACKED},
                  figs / "fig13_secant_history.png", title="secant iteration on the grillage")
    F.fig_mesh(list(MESH_NS), {f"w_pred, {n}": res["mesh"][n]["w_pred_mm"] for n in CRACKED}
               | {f"service deflection, {n}": res["mesh"][n]["deflection_service_mm"] for n in CRACKED},
               figs / "fig14_grillage_mesh.png", note="grillage mesh: relative change against N = 80")

    for name in CRACKED:
        r = res[name]
        print(name, "eta_dec", round(r["eta_decompression"], 4), "eta_ref", round(r["eta_ref"], 4),
              "table err", r["table"]["max_rel_error_at_midpoints"])
        for c in r["cases"]:
            print("  ", c["name"], c["status"], c.get("eta_mc") and {q: round(c["eta_mc"][q], 4) for q in QUANTILES},
                  c.get("eta_mc_exact_max_abs_diff"), c.get("eta_lower_bound_from_gate"))
    print("design", res["design"]["n_strands"], "timing", res["timing_s"])
    if pynite:
        print("pynite", res["pynite"]["max_abs_w_diff_mm"], res["pynite"]["max_w_mm"])


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_run_stage2.py tests/test_run_all_cases.py -q`
Expected: `5 passed` (2 + 3).

- [ ] **Step 5: Commit**

```bash
git add scripts/run_stage2.py tests/test_run_stage2.py
git commit -m "feat: run_stage2 - grillage inverse, Monte Carlo, indicator, PyNite, mesh, timing" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 12: Stage-2 results, progress log Entry 4, decisions, README, final verification

**Files:**
- Generate: `results/stage2.json`, `results/figures/fig9_grillage_w_eta.png`, `fig10_eta_bands.png`, `fig11_moment_share_G1.png`, `fig11_moment_share_G3.png`, `fig12_deflection_G1.png`, `fig12_deflection_G3.png`, `fig13_secant_history.png`, `fig14_grillage_mesh.png`
- Modify: `docs/PROGRESS_LOG.md`, `docs/DECISIONS.md`, `README.md`

- [ ] **Step 1: Run the Stage-2 script**

Run: `python scripts/run_stage2.py`
Expected (about 6 minutes): it prints, for G1 and G3, `eta_dec`, `eta_ref` and the table error (below 5e-3), one line per case (`measured` with the 5–50–95 % band and an exact-check difference below 1e-3; `refused` with a gate bound or `None` at GSD 0.38), the strand counts, the timing, and `pynite <max |Δw| in mm> <max |w| in mm>` with the first number below 1e-6 times the second. Probe values (without PyNite): G1 `eta_dec 0.514 eta_ref 0.2705`; G3 `eta_dec 0.6147 eta_ref 0.3528`; strands `(31, 29, 27, 29, 31)`. If the script stops on a non-converged secant loop (`build_table` raises "not finite"), do not change tolerances: rerun that analysis with `relax=0.5` in a scratch script to see whether it converges, and report the finding before changing the code.

- [ ] **Step 2: Check that only Stage-2 outputs changed**

Run: `git status --short results/` and `git diff --stat -- results/stage1.json results/figures/fig1_w_eta.png results/figures/fig2_mkappa.png results/figures/fig3_ei_ratio.png results/figures/fig4_deflection.png results/figures/fig5_eta_hist.png results/figures/fig6_mesh.png results/figures/fig7_cases.png results/figures/fig8_tornado.png`
Expected: the first lists only `results/stage2.json` and the eight new `fig9`…`fig14` files (untracked); the second prints nothing.

Look at every new figure (open the PNGs). Each must be readable: legends not covering data, labels in the plotted units, no empty panels.

- [ ] **Step 3: Append progress log Entry 4 and fill it from `results/stage2.json`**

Append this template to `docs/PROGRESS_LOG.md`, then replace every `<...>` with the value from `results/stage2.json` (keys in the table after the template) or from the test output. Write numbers as the log does (3–4 significant figures, units).

```markdown
## Entry 4 — B4–B5 Stage-2 grillage of S01 (<date>)

**Planned:** put the Stage-1 girders into a grillage of the whole S01 deck, so that the FE model sits inside the secant loop and the η inverse (spec P5, §3 Stage 2); η from a crack on the exterior girder G1 and on the interior girder G3, with a Monte Carlo band; load shed to the neighbours and the deflection of every girder (Plan 2).

**Found:**
- Grid element: reproduces Logan (5th ed.) Example 5.5 — reduced stiffness, d1y = −2.83 in, φ1x = 0.0295 rad, φ1z = −0.0169 rad and every printed local force — to Logan's three figures, and the sympy derivation exactly; the same example through assembly and partition gives the same numbers.
- Linear checks: identical girders under equal load each carry qL²/8 and deflect 5qL⁴/384EI; a lone girder line is the Stage-1 girder (one pass); a torsion-free deck follows Courbon's rule with rigid cross members and the lever rule with flexible ones, and the real deck lies between them. PyNite 2.0.2 agrees with the grid solve: max |Δw| = <dw> mm against a largest deflection of <wmax> mm; midspan moments to a relative <dM>.
- Design: strands G1–G5 = <n1>, <n2>, <n3>, <n4>, <n5> (Stage 1: 30) after <rounds> design rounds; soffit stress under service <s1> / <s2> / <s3> MPa (G1 / G2 / G3, all ≤ 0). Truck midspan moment at each girder's critical placement <Mt1> / <Mt2> / <Mt3> kN·m (G1 / G2 / G3); <which placement governs G1>.
- Blind zone: no crack at rest above η_dec = <edG1> (G1) and <edG3> (G3); Stage-1 interior girder <edS1>.
- η per measured case (Monte Carlo 5–50–95 %, 2,000 samples on the tabulated w(η); the exact inverse at the same widths differs by at most <dexact>):

| Case | η true (Stage-1 synthetic) | Stage 1 | Grillage, G1 cracked | Grillage, G3 cracked |
|---|---|---|---|---|
| case2_gsd0p02 | 0.36 | <band> | <band> | <band> |
| case2_gsd0p05 | 0.36 | <band> | <band> | <band> |
| case3_gsd0p02 | 0.21 | <band> | <band> | <band> |
| case3_gsd0p05 | 0.21 | <band> | <band> | <band> |
| case3_gsd0p10 | 0.21 | <band> | <band> | <band> |

- The synthetic widths were made with the Stage-1 interior girder, so "η true" is a Stage-1 truth: <one sentence on how far G1 and G3 move from it and why (fig9)>.
- Refusals: nominal gate bounds η ≥ <G1 bounds at GSD 0.02 / 0.05 / 0.10> (G1) and <G3 bounds> (G3); none at GSD 0.38.
- Indicator, G1 cracked at η = <etaG1> (reference case case2 at GSD 0.02): stage-B midspan moments G1–G5 intact <MiG1> → damaged <MdG1> kN·m; G1 loses <shed G1> % of its moment, G2 takes <shed G2> % and G3 <shed G3> % of it (fig11). Truck-only deflection of G1 <dti G1> → <dtd G1> mm (band <band G1>), limit L/800 = <lim> mm (fig12). Cracked elements per girder <crG1>; secant updates dead / service <itG1>.
- Indicator, G3 cracked at η = <etaG3>: moments <MiG3> → <MdG3> kN·m; G3 loses <shed G3>, G2 / G4 take <shed G2> / <shed G4> %; truck deflection of G3 <dti G3> → <dtd G3> mm (band <band G3>); cracked elements <crG3>; updates <itG3>.
- Secant loop on the grillage: more than one pass in every damaged analysis (P5 shown on the indeterminate model), monotone decrease of the moment change (fig13); <non-converged analyses, or "none">.
- Torsion sensitivity (GJ of cracked girder elements × 0.5): <G1: change in w_pred and in G1's midspan moment>; <G3: same>.
- Mesh (fig14), relative to N = 80: w_pred <G1 and G3 at N = 20>; damaged service deflection <G1 and G3 at N = 20>; the model uses N = 20.
- Table: largest relative interpolation error at the knot midpoints <tabG1> (G1), <tabG3> (G3).
- Timing: grillage build <t_build> s, linear solve <t_lin> ms, nonlinear service analysis <t_nl> s, table <t_table> s per girder, all cases with Monte Carlo <t_cases> s per girder; whole run <t_total> min.

**Limitations (for the report):**
- Every structural value is a stand-in (girder, strands, deck layout, overhang, diaphragms); the crack widths are synthetic and were generated with the Stage-1 girder.
- The diaphragms are their webs alone (not composite with the slab), a lower bound on their stiffness; bearings leave the girder ends free to twist.
- The initial pre-cracked zone comes from a linear intact envelope (the cracked girder's own shedding is ignored there, so it errs toward more pre-cracking); the truck stays at its Stage-1 longitudinal placement and at the critical transverse placement of the intact deck.
- The overhang's wet-slab torque on the bare exterior girder (stage A) is not modelled; stage A is vertical statics, as in Stage 1.
- No tornado on the grillage (DECISIONS 2026-10-15); the Stage-1 tornado stands.

**Changed & why:** <deviations from Plan 2 during execution, each with its DECISIONS line; "none" if none>.

**Numbers:** `results/stage2.json`; figures `results/figures/fig9…fig14`; `python -m pytest -q` → <summary line>.

**Reading guide:** read elements.py (grid part) → gridmodel.py → gridsolve.py → nonlinear.secant_iterate (the course part), then deck.py → grillage.py → wtable.py. The parts to defend in depth:

| Code | Course | Where |
|---|---|---|
| `elements.grid_k_local`, `grid_T`, `grid_k_global` | grid element: beam bending + torsion bar, transformation k = Tᵀ k̂ T | Logan (5th ed.) §5.4: Eq. 5.4.16 p. 241, Eqs. 5.4.17–5.4.18 p. 243, Eq. 5.4.19 p. 244; Example 5.5 pp. 244–250 |
| `gridsolve.stiffness`, `solve` | assembly member by member, partition K_ff D_f = F_f, reactions, local end forces f̂ = k̂ T d − f̂_eq | Ch3 slides 3-6 to 3-10; Logan §5.4 |
| `nonlinear.secant_iterate` | direct iteration with secant stiffness on an indeterminate model, D(i+1) = K_s(D(i))⁻¹ F | Ch11 slides, "Nonlinear Analysis Algorithms", 3) Secant Stiffness Algorithm |
| `deck.lever_weights`, `area_to_girders` | statically equivalent load transfer to the girder lines | statics (lever rule) |

Explain-back questions: (1) Why does the secant loop need more than one pass on the grillage when it needed exactly one on the single girder? (2) The same crack width gives a different η on G1, on G3 and on the Stage-1 girder (fig9, fig10): which model differences cause that? (3) When G1 cracks, G2 takes part of its moment: what does that do to G2's zero-tension design, and what would you check on the real bridge?

**Reviewed by student:** pending (after 2026-10-22).
```

Where the values come from:

| Placeholder | `results/stage2.json` key (or command) |
|---|---|
| `<date>` | the day Step 1 ran (YYYY-MM-DD) |
| `<n1>`…`<n5>`, `<rounds>` | `design.n_strands`, `len(design.design_rounds)` |
| `<s1>`, `<s2>`, `<s3>` | `design.soffit_stress_service_MPa[0..2]` |
| `<Mt1>`, `<Mt2>`, `<Mt3>` | `design.M_truck_mid_critical_kNm[0..2]` |
| `<edG1>`, `<edG3>`, `<edS1>` | `G1.eta_decompression`, `G3.eta_decompression`, `results/stage1.json` `eta_decompression` |
| case rows | `G1.cases[*]`, `G3.cases[*]` (`eta_mc.p05/p50/p95`, `eta_mc_exact_max_abs_diff`) and `stage1_comparison` |
| `<refusal bounds>` | `G1.cases[*].eta_lower_bound_from_gate`, `G3.cases[*].eta_lower_bound_from_gate` |
| `<etaG1>`, `<etaG3>` | `G1.eta_ref`, `G3.eta_ref` |
| moments, shares, shed | `G1.indicator.M_B_service_mid_kNm`, `.share_service`, `.shed_percent_of_Gk` (same for G3) |
| deflections | `G1.indicator.deflection_service_mm`, `.deflection_truck_mm`, `.deflection_limit_mm`; bands `G1.deflection_truck_band_mm` |
| cracked elements, iterations | `G1.indicator.cracked_elements`, `.iterations`, `.converged` (same for G3) |
| torsion sensitivity | `G1.gj_sensitivity`, `G3.gj_sensitivity` |
| PyNite | `pynite.max_abs_w_diff_mm`, `pynite.max_w_mm`, `pynite.max_rel_M_diff` |
| mesh | `mesh.G1`, `mesh.G3` (relative to N = 80) |
| table | `G1.table.max_rel_error_at_midpoints`, `G3.table...` |
| timing | `timing_s` |
| test counts | the summary line of Step 6 |

Run: `sed -n '/^## Entry 4/,$p' docs/PROGRESS_LOG.md | grep -n "<[^>]*>"`
Expected: no output (nothing in angle brackets remains in Entry 4).

- [ ] **Step 4: Append the execution decisions to `docs/DECISIONS.md` and fill them**

```markdown
## <date> — Plan 2 close-out (Task 12)

- **Stage-2 run:** reference η (case2 at GSD 0.02, Monte Carlo median) <etaG1> for G1 and <etaG3> for G3, against 0.361 on the Stage-1 girder; every analysis in run_stage2 converged (<secant update counts>) — or name the analysis that did not and what was done about it.
- **Neighbours under shed load:** cracked elements per girder <G1 case> (G1 cracked) and <G3 case> (G3 cracked); <one sentence: whether a neighbour cracks and what that means for its zero-tension design>.
- **Torsion kept uncracked on the main path:** GJ × 0.5 on cracked girder elements changes G1's w_pred by <%> and its stage-B midspan moment by <%> (G3: <%>, <%>).
- **Mesh n_el = 20 confirmed:** w_pred within <%> and damaged service deflection within <%> of N = 80.
- **PyNite agreement:** max |Δw| = <dw> mm, moments to <dM>; <differences found, or "none">.
- **Deviations from the plan's code during execution:** <each change with its reason, or "none">.
```

Run: `sed -n '/Plan 2 close-out/,$p' docs/DECISIONS.md | grep -n "<[^>]*>"`
Expected: no output.

- [ ] **Step 5: Update `README.md`**

Replace the paragraph that starts `Stage 1 (complete) models one interior girder` (it ends `grillage of the whole deck.`) with:

```markdown
Stage 1 (complete) models one interior girder of a simply supported span, built in two
construction stages: the bare girder carries the prestress and the weight of the girder and the
wet slab; the composite girder carries the wearing surface, the barriers and the design truck
through the beam finite-element model. Stage 2 (complete) puts the same girders into a grillage
of the whole five-girder deck (Logan's grid element, slab strips, diaphragms): the stiffness of
every girder element follows its own moment, so a cracked girder sheds load to its neighbours,
and the inverse for η runs through this indeterminate model for a crack on the exterior girder
G1 and on the interior girder G3.
```

In the Quickstart block, add after the line `    python scripts/run_all.py`:

```text
    python scripts/run_stage2.py
```

and replace the sentence that starts `` `run_all.py` regenerates `` with:

```markdown
`run_all.py` regenerates `results/stage1.json` and fig1–fig8 in `results/figures/`, and
`run_stage2.py` (about 6 minutes) regenerates `results/stage2.json` and fig9–fig14, both from the
committed pipeline outputs in `data/claimbridge/`; neither needs the private pipeline.
```

- [ ] **Step 6: Final verification**

Run: `python -m pytest -q`
Expected: `379 passed` with a pipeline configuration present (`376 passed, 3 skipped` without): the 308 Stage-1 tests and 71 new ones. Paste the summary line into Entry 4.

Run: `python -m pytest tests/test_stage1_golden.py -v`
Expected: `5 passed` — the Stage-1 fallback numbers did not move.

Run: `git status --short`
Expected: only `docs/PROGRESS_LOG.md`, `docs/DECISIONS.md`, `README.md` modified and the Stage-2 result files new. Run `git grep -n -i -E 'c:/users|c:\\users|/home/' -- . ':!*.png' ':!docs/superpowers/plans/*'` — expected: no output (no local paths committed; the plans are excluded because this command's own pattern is written in them).

- [ ] **Step 7: Commit**

```bash
git add results/stage2.json results/figures/fig9_grillage_w_eta.png results/figures/fig10_eta_bands.png results/figures/fig11_moment_share_G1.png results/figures/fig11_moment_share_G3.png results/figures/fig12_deflection_G1.png results/figures/fig12_deflection_G3.png results/figures/fig13_secant_history.png results/figures/fig14_grillage_mesh.png docs/PROGRESS_LOG.md docs/DECISIONS.md README.md
git commit -m "results: Stage-2 grillage run, progress Entry 4, decisions, README" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Do not push and do not merge to `main`: both are gates for the student (B6).

---

## Self-review (2026-10-15, at plan writing)

Spec coverage (spec §3 Stage 2, §5.4, §5.5, §6 grillage rows, controller list of 2026-10-15):

| Requirement | Task |
|---|---|
| Grid element, 3 DOF/node, bending = `beam_k`, torsion GJ/L, [T], consistent loads; sympy + Logan example | 2 (test_grid_element), 3 (example through assembly) |
| Deck: 5 girders at 2.2 m, strips at every girder station, end + midspan diaphragms, overhang; Hambly properties; ν = 0.2 | 1 (inputs, S22–S27), 5 |
| Stage A per girder by statics with its own tributary width; stage B on the grillage (wearing by lever rule, barriers on the exterior girders with torque, truck as wheel loads with (1 + I) × 1.30, transverse scan with the 2 ft curb rule, longitudinal max-midspan placement, no distribution factor) | 5, 6 |
| Nonlinear: per-element staged sections and EC2 secant EI from the own stage-A moment; others η = 1 under the same rule; transverse/torsion linear; GJ × 0.5 sensitivity; sticky flag; stiffness residual; under-relaxation; iteration cap; non-converged result object | 4, 6, 7, 11 |
| Strand design per girder by the zero-tension rule under its own grillage service moment, tested | 6 |
| Inverse + Monte Carlo through a table, exact confirmation of p05/p50/p95; data reused unchanged | 9, 11 |
| Outputs: stage2.json; η bands G1/G3 vs Stage 1; moment share; deflections; load shed; convergence history; mesh convergence; PyNite | 8, 10, 11, 12 |
| Verification bar: element vs sympy and Logan; one girder = Girder1D; identical girders, ΣR, K symmetric; lever rule / Courbon; truck share vs S/5.5; PyNite; >1 pass with consistent EI; inverse round trip | 2, 3, 5, 6, 7, 8 |
| Girder1D and Stage 1 untouched, golden test guards | every task (append-only to Stage-1 files; golden test in Tasks 1, 4, 6, 12) |
| Repo hygiene; no push | conventions, Task 12 Step 6 |
| Tornado on the grillage (spec §8 B5) | not in this plan — see DECISIONS (Task 1) |

Placeholder scan: the only angle-bracket placeholders are inside the Entry 4 and close-out templates of Task 12, which the controller asked for and which Step 3/Step 4 check are all filled. Name consistency checked against the real code: `Girder1D(inp, n_el, dy_max, n_strands=...)`, `.sec_el`, `.sec_nd`, `.M_A_el`, `.M_A_nd`, `.mcr(eta, sec, fct)`, `.crack_at(eta, M, sec)`, `.eps_p0(eta)`, `.EI_I`, `.make_section`, `.is_girder`, `.h_girder`, `.stiffness_fn(eta, precracked)`, `.precracked_zone(eta)`, `.analyse(eta, model)`, `.model_service`, `.q_b`; `design_prestress_staged(make_section, is_girder, h_girder, strand_area, fpe, M_A, M_total)`, `prestress_staged(..., n_strands, strand_area, fpe, M_A, M_total).soffit_stress_service`; `ec2_secant_EI(sec, eps_p0, M, Mcr, beta, M_ref)`; `estimate_eta(w_pred, w_meas)` → `EtaEstimate(status, eta, eta_lower_bound, iterations, evaluations)`; `eta_decompression(w_pred, eta_lo, eta_hi, w_tiny, xtol)`; `monte_carlo_eta(w_pred, w_meas, u95, n)` → `.summary()`; `truck_position_max_midspan(inp)` → (positions, moment); `axle_offsets(inp, reverse)`; `impact_factor(inp)`; `truck_fraction(inp)`; `governing_width(parsed)`, `refused_below_resolution(parsed)`; `figures.fig_mesh(Ns, series, path, note)`. Every Python block was compiled (`compile()` after dedent) when this plan was written.
