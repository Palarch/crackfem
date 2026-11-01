# Design — Crack-Width-Informed Prestress & Load-Sharing Indicator (CPK72 S01)

- **Course:** 01203541 Finite Element Analysis of Structures, KU, 1st sem. 2026 — solo project (15%, CLO9), constraint "AI + FEM"
- **Student:** Palid Archavachalee (6914500959)
- **Date:** 2026-10-07 · **Due:** 2026-11-01 · **Hand-in:** GitHub repo link + report (docx + PDF, English)
- **Status:** design approved in chat 2026-10-07; implementation plan to follow
- **Baseline:** the approved proposal "FEM–AI Structural Indicator for Crack-Width-Informed Deflection: A Case Study on Bridge CPK72 (Route 304)", dated 15 Sep 2026, approved without comments

---

## As-built amendments (2026-10-12)

The approved text below is not rewritten. Where the as-built Stage 1 differs from it, this list governs. Every item is a dated line in `docs/DECISIONS.md` (heading named in brackets).

1. **Staged construction adopted** (student decision) — §5.3's "no construction staging" is superseded: stage A (prestress, girder, wet slab) acts on the bare girder, stage B (wearing surface, barriers, truck) on the composite girder. The secant EI of §5.2 item 4 is measured from the stage-A moment, not from the prestressed, unloaded state. [DECISIONS 2026-10-12 — Student decisions; G09 review (the size of the staging shift)]
2. **Live load HS20-44 × 1.30** (student decision, JICA 2000 "DOH standard") — §5.3 names plain HS20-44. [2026-10-12 — Student decisions]
3. **Stage 1 uses interior-girder rules** (tributary width S, distribution S/5.5); the exterior girder G1 of §5.4 is handled in Plan 2. [2026-10-12 — Student decisions]
4. **Held-out cracked-zone check deferred to field data** (§5.4): the synthetic masks are drawn from the same model, so the comparison would be circular; the predicted zone is reported as a model output. [2026-10-12 — Plan 1 close-out (Task 22)]
5. **Tornado parameters as built** (§5.5): girder spacing, strand centroid, slab thickness, bond factor ξ, strand cover, girder fck, with the as-built strand count held (not re-designed per input). [2026-10-12 — Student decisions (ξ); G09 review (strand count held)]
6. **Where concrete physics meets FE** (§4): `girder1d.stiffness_fn` builds the per-element secant EI from the section models and hands it to `nonlinear.solve_secant`, which stays generic (it knows only EI as a function of the element moments). [2026-10-12 — Stage-1 close-out (final review)]
7. **Bridge to the thesis pipeline** (§4): the pipeline CLI prints a text report; `claimbridge.py` parses it into a dict and the saved files are JSON (parsed values plus the recognised report lines). [2026-10-12 — G09 review (parser fails closed); Stage-1 close-out]
8. **File names as built** (§4, §7): `scripts/make_synthetic_cases.py` (spec: `make_masks.py`), `src/crackfem/claimbridge.py` (spec: `io_claimbridge.py`), `scripts/run_claimbridge.py` (spec: `run_claimbridge.md`, a documented manual step). [2026-10-12 — Stage-1 close-out (final review)]
9. **PyNiteFEA** is approved (§7, §10) but not yet installed; `requirements.txt` pins the six packages Stage 1 uses, and PyNiteFEA is pinned when Plan 2 needs it. [2026-10-12 — Stage-1 close-out (final review)]
10. **Branson cross-check uses the same EC2 section and Mcr** (§5.1, §5.2 item 3, §4): the Branson interpolation replaces only the EC2 §7.4.3 interpolation; the ACI formulas for Ec and fr are not used in any analysis, and Bischoff's form is not implemented. [2026-10-12 — Stage-1 close-out (final review)]

## As-built amendments — Stage 2 (2026-10-15)

Same rule as above: the approved text is not rewritten; where the as-built Stage 2 differs, this list governs, and every item is a dated line in `docs/DECISIONS.md` (heading named in brackets; "design" = 2026-10-15 — Plan 2 design, "close-out" = 2026-10-15 — Plan 2 close-out (Task 12)).

1. **Cracked girders in the grillage** (§5.4): G1 (exterior) and G3 (interior), each analysed with the other four girders at η = 1 under the same secant rule. [design — Cracked girders; Nonlinear grillage]
2. **PyNiteFEA 2.0.2 pinned** (§7, §10; supersedes Stage-1 item 9): the newest release that keeps the pinned numpy; it cross-checks the linear intact grillage only. [design — PyNiteFEA pinned at 2.0.2; close-out — PyNite agreement]
3. **Modules as built** (§4: "`model.py` … grillage"): the grillage lives in `deck.py` (layout, Hambly properties, load mapping, truck placements), `gridmodel.py` (model data), `gridsolve.py` (assembly and partition solve, end forces) and `grillage.py` (per-girder sections, strand design, nonlinear analysis, w_pred, indicator); plus `wtable.py`, `pynite_check.py` and `scripts/run_stage2.py`. `model.py` and `solve.py` stay the 1-D girder. [design — Grid element; Deck model]
4. **Stage A by statics inside the grillage** (§5.3): prestress, girder and wet slab act on each bare girder (its own `Girder1D`); only stage B (wearing surface, barriers, truck) goes through the grillage. Stage A therefore stays on the cracked girder whatever its η: prestress lost after composite action is not restrained by the deck (this pushes η lower, the same direction as the shedding). [design — Stage A per girder by statics; close-out]
5. **Pre-crack flags and cracking path** (§5.2): the initial pre-cracked zone comes from a linear intact moving-truck envelope; stage B is applied incrementally and new cracks are flagged one event at a time (largest M/Mcr first) from converged states only, so the flags do not depend on the iteration path or the load-step size. [design — Initial pre-cracked zone; Incremental stage-B loading; Flags set one event at a time]
6. **Monte Carlo on a tabulated w(η)** (§5.5, §2 "no surrogate"): w_pred on a 99-point η grid plus the plateau edge, joined linearly — an interpolation of the exact map, not a trained surrogate; the point estimate is exact and the p05/p50/p95 are confirmed by the exact inverse (≤ 1.4e-4 in η). [design — Inverse and Monte Carlo through a table]
7. **Label `eta_gen_1d`** (§5.4, §6): the synthetic widths were generated with the Stage-1 single-girder model, so the η they were made with is reported as "η used to generate the crack (1-D model)", never "η true", in every Stage-2 output. [close-out — Labelling]
8. **Truck placement and diaphragms** (§5.3, §5.4): two trucks in non-overlapping 12 ft AASHTO design lanes (default; "clearance" as a sensitivity); diaphragm model "web" (default), "none" or "composite" as a sensitivity. [design — Lane rule; Diaphragm three-case sensitivity; close-out — G1 model sensitivity]
9. **No tornado on the grillage** (§8 B5): the Stage-1 tornado stands; the grillage has the G1 model sensitivity and the GJ × 0.5 run instead. [design — Outputs]
10. **Load shed defined** (§6): `shed_percent_of_Gk[g]` = (M_damaged,g − M_intact,g) / M_intact,k × 100, with M the stage-B service midspan moment (dead + truck at Gk's critical placement) and k the cracked girder; the entry for g = k is minus the share Gk loses ("G1 loses 71.6 %" = −shed[k]). Stage A is not included (it is unchanged). [close-out; 2026-10-15 — Stage-2 final review]

---

## 1. Goal

Turn a crack width measured by the author's thesis inspection pipeline (ClaimBridge, private) into a **mechanical** statement about the girder and the deck, which a categorical condition rating cannot give:

> "Girder Gk shows a w mm flexural crack at rest → it has about η of its prestress left (5–95% band) → under the design truck its midspan deflection is δ, its stiffness ratio along the span is Ie/Ig(x), and its neighbours now carry X% more moment."

The result is an **advisory structural indicator** shown next to the rule-engine rating for the same crack. It never feeds back into that rating.

FEM must be the engine, not a post-processor: the girder's stiffness depends on its own moment, and in an indeterminate deck that moment depends on every girder's stiffness. So the answer cannot be obtained without solving the FE model inside a nonlinear loop.

## 2. What changed from the proposal, and why

The approved proposal is the baseline. Five issues found during design (2026-10-06/07) drive the changes. Each one becomes a chapter of the progress report.

| # | Proposal said | Problem | Change |
|---|---|---|---|
| P1 | Crack width w sets Ie through Branson's equation | Branson has no w term; it uses Ma/Mcr. The w → stiffness link was a placeholder curve ψ(w) | Explicit chain: w → strand stress increase Δσp (inverted EC2 §7.3.4) → section state → stiffness |
| P2 | Damage in one element; mesh convergence N = 10/20/40 against closed form | Damaged length = L/N, so the result depends on the mesh (2.45 vs 0.61 m). Also, Hermite beam elements with consistent loads give exact nodal deflections for any N, so matching the closed form proves nothing | Stiffness assigned per element from that element's own moment; the cracked zone emerges where M > Mcr. Convergence measured on the continuous EI(x) |
| P3 | Equivalent RC rectangle 0.45 × 1.70 m, RC Branson | A 24.5 m highway span is typically a precast prestressed I-girder (as-built drawings not available) | Stand-in standard PC I-girder + composite slab, prestress sized by a design rule |
| P4 | Crack width "comes straight from ClaimBridge" | The pipeline has not yet measured a real crack (first field survey: Nov 2026) | Synthetic crack masks of known width run through the real crack-width command; stated plainly in the report |
| P5 | Iterate Ie ↔ Ma until deflection stops changing | S01 is simply supported, so M(x) comes from statics and never changes: the loop ends after one pass | Model the deck as a grillage (indeterminate): the loop becomes a real nonlinear solve. The 1-D stage keeps a test proving the one-pass convergence |

Also dropped: the proposal's "future work" ML surrogate. A beam/grillage solve costs milliseconds, so a surrogate would add nothing. The report shows the measured timing as the evidence.

## 3. Scope and stages

- **Stage 1 — 1-D girder (fallback deliverable).** Whole chain on one girder: section model, crack width, beam FE, nonlinear loop (provably 1 pass), inverse for η, Monte Carlo, ClaimBridge link, figures. Complete and submittable on its own.
- **Stage 2 — grillage of S01 (main deliverable).** Same chain, with `model.py` producing a multi-girder grillage. The nonlinear loop now redistributes moment between girders.

**Out of scope:** real CPK72 section/strands; full moving-load envelopes (only the critical truck position(s)); creep, shrinkage and relaxation beyond what η absorbs; shear cracking; dynamics; 2-D continuum models; any change to ClaimBridge code; an ML surrogate (see §2).

## 4. Architecture

```
[1] scripts/make_masks.py   synthetic crack masks (known width, known cracked-zone length)
        │ PNG mask + surface JSON (chosen GSD)
        ▼
[2] ClaimBridge crack-width command (unchanged, run locally)
        │ w ± u (mm) and rule-engine rating → saved to data/claimbridge/*.json
        ▼
[3] crackfem (Python package)
    params.py       every numeric input, each with a source tag (see docs/sources.md)
    section.py      elastic fibre section of the composite PC girder → curvature, strand stress, Mcr(η)
    crackwidth.py   EC2 §7.3.4 crack width for bonded prestressed members; closed-form inverse w → Δσp
    stiffness.py    EC2 §7.4.3 interpolation → secant EI per element; ACI Branson/Bischoff as a cross-check option
    elements.py     beam element (Ch3) and grid element (Logan §5.4): k, [T], consistent loads
    model.py        S01 geometry → nodes, elements, DOF map, supports (1-D girder | grillage)
    solve.py        assembly + partition solve (course method); member end forces
    nonlinear.py    inner loop: EI_e = f(M_e, η) → secant/direct iteration until moments converge
    inverse.py      outer loop: safeguarded secant/Newton on η until w_pred(η) = w_meas
    uq.py           Monte Carlo on w; one-at-a-time sensitivity (tornado)
    io_claimbridge.py  reads the pipeline's output JSON
    figures.py      all figures from results JSON
        ▼
[4] results/*.json + results/figures/* → report
```

**Interface rules**
- `elements`, `solve` and `model` know nothing about concrete or cracks. They take EI (and GJ) per element and return displacements and member forces.
- `nonlinear.py` is the only module where concrete physics meets the FE model.
- ClaimBridge and Python exchange only JSON files. No ClaimBridge code or data other than the synthetic masks and their outputs enters this repo.

## 5. Physics and stand-in assumptions

All values below are **candidates**. Block B1 (source pack) opens each source, records clause and page in `docs/sources.md`, and locks the value in `params.py`. Paywalled originals (EN 1992-1-1, ACI 318-19, Hambly) are replaced by open secondary sources that quote the clause, tagged `[secondary]`.

### 5.1 Stand-in structure (every item labelled "stand-in" in outputs)

| Item | Candidate | Source to verify |
|---|---|---|
| Span | S01, L = 24.50 m, simply supported | proposal; drawing reference already used by ClaimBridge |
| Girder | AASHTO Type IV PC I-girder (≈ 1.37 m deep) | PCI Bridge Design Manual property table |
| Deck | 5 girders at ≈ 2.2 m, 0.20 m composite slab, end diaphragms + one midspan diaphragm | none — parameter; covered by sensitivity |
| Strands | 12.7 mm seven-wire low-relaxation, fpu = 1860 MPa | ASTM A416 / PCI |
| Prestress amount | Sized so the intact girder (η = 1) has zero bottom-fibre tension under full service load | AASHTO Class U / zero-tension design rule |
| Concrete | girder fck 40 MPa, slab 28 MPa; Ecm and fctm from EC2 Table 3.1 (main path); ACI Ec and fr only inside the Branson cross-check | EN 1992-1-1 Table 3.1; ACI 318-19 §19.2 |
| Grillage rules | longitudinal members = composite girders; transverse members = slab strips (I, J per unit width) + diaphragms | Hambly, *Bridge Deck Behaviour* |
| Units | SI throughout (kN, m, MPa) | — |

### 5.2 Physics chain

1. **Section** — elastic fibre section: concrete linear in compression and zero in tension once cracked; strands bonded and linear elastic with initial prestrain from η·P_design. Returns curvature and Δσp for a given moment; Mcr(η) comes from the same model and is checked against P/A + Pe/S + fr.
2. **Crack width** — EC2 §7.3.4: wk = sr,max × (εsm − εcm), prestressed version with Δσp and bond factor ξ1. Linear in Δσp, so w → Δσp is closed form (including the 0.6·Δσp/Es floor branch).
3. **Stiffness** — EC2 §7.4.3: curvature = ζ·κ_cracked + (1 − ζ)·κ_uncracked, ζ = 1 − β(Mcr/M)², β = 0.5 (sustained or repeated loading, S10 — a bridge girder sees both; changed from 1.0 on 2026-10-07). Same code family as the crack width. ACI Branson/Bischoff run as a cross-check.
4. **FE convention** — FE carries external loads only, measured from the prestressed, unloaded state: secant EI_e = M_e / (κ(M_e) − κ(0)). Camber is reported separately.

### 5.3 Load states

- **Inspection state (fits η):** dead load only — girder, slab, 50 mm wearing surface, barriers. Simplification: all dead load acts on the composite system (no construction staging); stated in the report.
- **Indicator state (outputs):** dead load + HS20-44 truck positioned for maximum moment in Gk (verify the DOH design-load reference in B1).

### 5.4 Inverse problem

- Unknown: η ∈ (0, 1], the remaining prestress fraction in the cracked girder Gk (other girders intact).
- Cracked girder: default Gk = G1, the exterior girder (most exposed to deck drainage — an inference, labelled as such); interior G3 is run as a second case. In Stage 1 the single girder is Gk.
- Observation fitted: the largest per-crack maximum width that the pipeline reports for Gk's mask. The model evaluates it at Gk's maximum-moment section (the cracks in the mask are drawn centred on midspan).
- Forward map w_pred(η): nonlinear grillage solve (inspection state) → M at the critical section of Gk → section → Δσp → EC2 crack width.
- Solver: safeguarded secant on a bracket (falls back to bisection). Monotonicity of w_pred(η) is checked on a grid first; if it is not monotonic, the report says the problem is ill-posed for that case.
- w below the cracking threshold → returns the bound "η ≥ η_decompression", never a point value.
- **Held-out check (not fitted):** predicted cracked-zone length (where M > Mcr(η) along Gk) vs the crack extent drawn in the synthetic mask.

### 5.5 Uncertainty

- Monte Carlo on w only: w ~ Normal(w_meas, (u/2)²), truncated at w > 0, with ClaimBridge's ±u read as a 95% bound. Fixed seed; n = 2,000 by default, with a stability check at 1k vs 10k on the 1-D stage.
- One-at-a-time sensitivity (tornado) on the main stand-ins: girder spacing, f'c, cover, sr,max, slab thickness.
- Outputs: η, δ_mid(Gk), Ie/Ig(x), moment share per girder vs intact — each as median and 5–95% band.

### 5.6 Known limits (go into the report as stated)

- A crack seen at rest may be a residual opening from a past overload (debris prevents closure) → η would be underestimated.
- Every structural value is a stand-in until the BRB-CP drawings are available.
- Crack widths are synthetic until the Nov 2026 field survey. The thesis pipeline is described only as far as its code supports.

## 6. Verification

Every check compares against evidence independent of the code under test.

| Module | Check | Independent evidence |
|---|---|---|
| elements | beam k and consistent loads | symbolic derivation in sympy from the shape functions (k = ∫ Bᵀ EI B dx) |
| elements | grid element (torsion + bending + [T]) | reproduce a Logan §5.4 worked example to the printed answers |
| solve | SS beam, uniform load, N = 1…40 | 5wL⁴/384EI exact at nodes for every N (recorded as a finding, see P2) |
| solve | stepped EI beam | virtual-work integral with scipy quadrature |
| grillage | 1 girder, no cross members | equals the 1-D beam result |
| grillage | identical girders, uniform load | equal shares; ΣR = total load; K symmetric; response symmetric |
| grillage | truck load share | within a sanity band of the AASHTO distribution factor |
| grillage | independent solver | same model in PyNite (pinned version) |
| section | Mcr(η) | hand calculation P/A + Pe/S + fr |
| section | cracked prestressed section | a published worked example (located in B1) |
| crackwidth | EC2 §7.3.4 | a published EC2 worked example; w → Δσp → w round trip incl. floor branch |
| nonlinear | 1-D girder | converges in exactly one pass (proves P5) |
| nonlinear | grillage | moment change < 0.1%; residual history saved; at η = 1 equals the linear solve |
| inverse | round trip | η_true → w → recovered η for several values; monotonicity grid |
| inverse | no-crack case | returns the bound, not a number |
| ClaimBridge link | synthetic mask of known width | the pipeline's result within ClaimBridge's own accuracy gate; its own crack-width test suite run locally first |
| uq | zero uncertainty | collapses to the point answer; percentiles stable 1k vs 10k |
| uq | timing | measured solve time (evidence for dropping the surrogate) |

Report level: every figure is drawn by a script from saved results JSON; every number in the text traces to that JSON; the final `pytest` output goes in the appendix.

## 7. Repository layout

```
fem_project/
  README.md                 3-command quickstart
  requirements.txt          pinned (numpy, scipy, sympy, matplotlib, PyNiteFEA, pytest, pillow)
  src/crackfem/             modules from §4
  tests/                    one test file per module + end-to-end
  scripts/                  make_masks.py, run_claimbridge.md (documented manual step), run_all.py
  data/claimbridge/         synthetic masks + the pipeline's output JSON (committed)
  results/                  generated JSON + figures (regenerated by run_all.py)
  docs/
    PROGRESS_LOG.md         one entry per block: planned → found → changed & why → numbers → reading guide → what the student reviewed
    DECISIONS.md            one dated line per non-obvious decision
    sources.md              every source: clause/page, open or [secondary], value used
    superpowers/specs, superpowers/plans
  report/                   docx + PDF (B6)
```

## 8. Timeline

| Dates | Block | Student time |
|---|---|---|
| Oct 7–8 | spec + implementation plan | read spec |
| Oct 8–10 | **B1** source pack → `sources.md`, `params.py` | gate answers only |
| Oct 10–13 | **B2** Stage 1 core (section, crack width, stiffness, beam FE, nonlinear, inverse) + tests | 0 |
| Oct 13–16 | **B3** ClaimBridge link + Monte Carlo on 1-D + first figures → fallback complete | 0 |
| Oct 16–19 | freeze (mock + final exam) | 0 |
| Oct 20–23 | **B4** Stage 2 grillage: grid element, deck model, nonlinear loop, PyNite cross-check | ~1 h/day review |
| Oct 24–26 | **B5** inverse + Monte Carlo on grillage, tornado, final figures | ~1 h/day review |
| Oct 27–30 | **B6** report from PROGRESS_LOG, review sessions, fixes | ~1.5 h/day |
| Oct 31 | buffer | |
| Nov 1 | submit | |

Slip rule: if B4 is not done by Oct 25, submit Stage 1 complete + the grillage as a section marked "partial".

## 9. Authorship and review

The code is written by Claude (Anthropic's AI) and reviewed by the student. Each block's log entry ends with a reading guide (file order, what each function does in course terms, three explain-back questions) and marks the parts the student must be able to defend in depth: the beam/grid element, assembly/partition, and the secant loop. The report states the authorship model plainly and includes the student's review notes.

## 10. Gates and risks

- **Dependencies:** matplotlib and PyNiteFEA, pinned — approved 2026-10-07.
- **GitHub push:** asked at the time, including public vs private. The repo must contain nothing private from the thesis pipeline (no field data, no internal URLs or keys) — only synthetic masks and their outputs.
- **Nothing published, nothing written into the ClaimBridge repo.**
- **Risks:** secondary sources may not quote every clause needed (then the report states the gap); the pipeline may require real surface geometry instead of a synthetic surface JSON (B3 decides; real field data is never committed); grillage schedule is tight after the exam (slip rule above).
