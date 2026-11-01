# Progress Log — crackfem

Each entry: **Planned → Found → Changed & why → Numbers → Reading guide → Reviewed by student.**
Numbers come only from `results/*.json` or test output. The code is written by Claude (Anthropic) and reviewed by the student; the "Reviewed by student" line records what was reviewed and when.

## Entry 0 — Design (2026-10-06 → 2026-10-07)

**Planned:** execute the approved proposal (1-D Euler–Bernoulli beam of CPK72 S01, Branson Ie at one cracked element, ±2·GSD band from ClaimBridge).

**Found:** five problems in the proposal (spec §2):
- P1: Branson's equation has no crack-width term — the w → stiffness link was a placeholder curve.
- P2: damage forced into one element makes the result mesh-dependent, and for the uncracked beam Hermite elements are exact at the nodes for any N, so the planned convergence check could not prove anything.
- P3: a 24.5 m highway span is typically a prestressed I-girder, not RC.
- P4: the thesis pipeline has not yet measured a real crack (field survey: Nov 2026).
- P5: a simply supported span is statically determinate, so the planned Ie ↔ Ma iteration ends after one pass.

**Changed & why:** inverse for remaining prestress η from the crack width; EC2 crack width + stiffness; per-element stiffness from each element's own moment; 1-D girder (Stage 1) → grillage (Stage 2) so FE sits inside the nonlinear and inverse loops (spec §2–§5, DECISIONS 2026-10-07).

**Numbers:** none yet. A probe on 2026-10-07 confirmed the pipeline's crack-width command measures a synthetic 6 px crack at GSD 0.05 mm/px as 0.30 ± 0.10 mm and the rule engine rates it 5 (Excellent).

**Reading guide:** spec §1–§5.

**Reviewed by student:** spec approved in chat 2026-10-07.

## Entry 1 — B1 source pack (2026-10-07)

**Planned:** open every source behind params.py (spec §5) and lock the values.
**Found:**
- S01 span 24.50 m: [primary], BRB-CP sheet 3.
- S02 deck layout: [stand-in]; barrier weight now from TxDOT T501 (326 plf/side); Thai I-girder spacing of about 2.25 m (Kasem Bundit Eng. J. 2016) supports 2.20 m.
- S03 Type IV: [primary] PCI BDM dimensions; three published I values (PCI 260,730; NCDOT 260,741; TxDOT 260,403 in⁴); NCDOT equals the exact outline integral.
- S04–S08, S10, S13, S14, S16–S18, S20: [secondary] only. The paywalled EN 1992-1-1, ACI 318 and AASHTO Standard Specifications were read through JRC, European Concrete Platform, Sétra, TxDOT/TTI, PCI Journal, NCHRP and software-vendor documents that quote them.
- S09 ξ = 0.6 (pretensioned strand, Table 6.2): unverified; only Eq. 7.5 (ξ1 = √ξ) and the post-tensioned row (0.5) were found in open sources.
- S11 unit weights: [secondary]; AASHTO 150/140 pcf vs EN 1991-1-1 25 kN/m³ vs JRA 24.5 kN/m³ disagree; AASHTO adopted.
- S12 HS20-44: [secondary]; a JICA (2000) design table gives the DOH standard as HS20-44 × 1.30, not applied.
- S15 strand centroid/cover and S19 concrete strengths: [stand-in], with reasoning from the PCI 2-in strand grid and the TxDOT Type IV example.
- S21: EC2 crack-width example found (ECP Ex. 7.3, wk = 0.184 mm). Cracked prestressed example: Mast 1998 Ex. 1, entered with the printed decompression stress fdc (the model's locked-in strain is eps_p0 = fdc/Ep).
**Changed & why:**
- barrier_weight 7.0 → 4.76 kN/m (S02).
- unit_weight_concrete 24.0 → 23.56 kN/m³ (S11).
- unit_weight_wearing 22.0 → 21.99 kN/m³ (S11, now stated as 140 pcf).
- All three are now traceable to an opened source, and the dead load is in the same AASHTO family as the live load (DECISIONS 2026-10-07, B1).
**Numbers:** see docs/sources.md.
**Reading guide:** skim sources.md; rows tagged [stand-in] are the assumptions the examiner may question.
**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 2 — B2 Stage-1 core (2026-10-07 → 2026-10-12)

**Planned:** section, crack width, stiffness, beam FE, nonlinear loop, inverse — each with independent tests (spec §6).
**Found:**
- Beam nodal deflections are exact for N = 1, 2, 3, 5, 8 and 40 elements (proposal P2 confirmed: a mesh study on the uncracked beam proves nothing).
- The secant loop stops after one pass on a determinate beam (P5 confirmed), including the staged Stage-1 girder with damage (η = decompression edge − 0.15, iterations_service = 1).
- An indeterminate toy beam (propped cantilever, softening stiffness) needs at least 2 passes and sheds fixed-end moment to the span; there the loop is real.
- A section whose strand is below the kern cracks at the top under prestress alone, because cracked mode lets concrete carry compression only anywhere (test corrected in Task 7; curvature references come from uncracked states).
- The staged section matches the two-stage transformed-section hand calculation (soffit stress, stage-B curvature and composite EI within 1e-5), and a composite solve at the stage-A moment leaves the slab unstressed (below 1 Pa).
- Published examples reproduced: EC2 crack width (ECP Ex. 7.3) and the cracked prestressed section (Mast 1998 Ex. 1).
- Stage-1 girder: the intact girder shows no crack at rest, w_pred does not increase with η on a 40-point grid, the inverse recovers η through the real section to 1e-5, and the intact service deflection equals the linear uncracked solution.
**Changed & why:** staged construction adopted (student decision 2026-10-12; measured shift +0.06 in the soffit-decompression η, +0.10 for the same girder; an earlier estimate of 0.15 mixed two inconsistent assumptions — see DECISIONS 2026-10-12 G09 review); HS20-44 × 1.30 applied (student decision 2026-10-12); the truck scan adds the exact kink placements (an axle at midspan or at a support) because the 1 mm grid can miss the peak; the secant loop also checks the stiffness residual before stopping, because a moment-only criterion can stop with an EI that belongs to other moments when stiffness jumps at Mcr; the inverse finds the decompression edge by bisection, checks the residual after Illinois so a jump in w_pred gives a bound instead of a fake point estimate, treats w ≤ 1 nm as no crack and refuses NaN; two plan test errors corrected in Task 13 (round-off bound 1e-9 instead of 1e-12, vertical reactions only in the force check) — all in DECISIONS 2026-10-12.
**Numbers:** `python -m pytest -q` → 165 passed, 0 skipped.
**Reading guide:** read in this order: elements.py → model.py → solve.py → nonlinear.py (the course part), then section.py → staging.py → crackwidth.py → stiffness.py → girder1d.py (the physics). Explain-back: (1) why are the nodal deflections exact for any N? (2) why does the secant loop stop after one pass here? (3) why does the slab carry no stress from its own weight, and what does that change in the crack width at rest?
**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 3 — B3 thesis-pipeline link, Monte Carlo, Stage-1 results (2026-10-12)

**Planned:** run synthetic cracks through the real thesis pipeline, invert for η with a Monte Carlo band, verify mesh convergence and the Branson cross-check (spec §5.4–§6).

**Found:**
- Design (staged; smallest strand count with no soffit tension under service at η = 1): 30 strands, P0 = 3212 kN (fpe = 1085 MPa).
- Blind zone: no crack is predicted at rest for η ≥ η_dec = 0.557, so a loss of up to 44 % of the prestress leaves no visible sign. Below η_dec the predicted width does not increase with η and reaches 0.89 mm at η = 0.02 (fig1).
- Resolution: the thesis pipeline refuses cracks narrower than 3 px. At GSD 0.38 mm/px (the lower end of its field GSD) every case was refused as "below-measurable-resolution": the gate (3 px = 1.14 mm) is wider than any crack this girder can show (0.89 mm), so nothing is measurable and a refusal gives no bound. In all, 7 of 12 cases were refused (fig7). A refusal at a finer GSD still bounds η from below through the 3 px width: nominal η ≥ 0.439 / 0.368 / 0.272 at GSD 0.02 / 0.05 / 0.10 mm/px; one-sided 95 % (σ = 1 px) η ≥ 0.411 / 0.313 / 0.179.
- Inverse accuracy, the five measured cases (2,000 Monte Carlo samples each, none dropped):

| Case | GSD (mm/px) | Mask width (px) | w ± u (mm) | η 5–50–95 % | η true | DOH rating |
|---|---|---|---|---|---|---|
| case2 | 0.02 | 8 | 0.16 ± 0.04 | 0.338–0.361–0.385 | 0.36 | 5 |
| case2 | 0.05 | 3 | 0.15 ± 0.10 | 0.312–0.367–0.433 | 0.36 | 5 |
| case3 | 0.02 | 20 | 0.40 ± 0.04 | 0.196–0.215–0.234 | 0.21 | 5 |
| case3 | 0.05 | 8 | 0.40 ± 0.10 | 0.169–0.214–0.262 | 0.21 | 5 |
| case3 | 0.10 | 4 | 0.40 ± 0.20 | 0.124–0.214–0.312 | 0.21 | 5 |

- η true lies inside every band, and the band widens with GSD because u = 2·GSD (case3: 0.037 → 0.093 → 0.188 wide). The bands are measurement noise only, conditional on the model (fig1, fig5).
- Rule engine: every measured case is rated 5 "ดีมาก (Excellent)" — 0.16 mm, where the model says η ≈ 0.36 (64 % of the prestress lost), and 0.40 mm, where it says η ≈ 0.21 (79 % lost). The rating carries no information about prestress; the indicator sits next to it and never changes it.
- Serviceability at η_ref = 0.36 (case2 at GSD 0.02, the first detectable case): truck-only midspan deflection 72.1 mm against 9.5 mm intact (fig4); Monte Carlo band 69.1–72.1–74.4 mm (5–50–95 %, read from η p95–p50–p05 in reversed order, because deflection falls as η rises). Limit L/800 = 30.6 mm: the intact girder uses 31 % of it; the damaged girder deflects 2.4 times the limit. Secant EI / uncracked EI falls to 0.12 at midspan; 24 of 40 elements (x ≈ 4.9–19.6 m) are pre-cracked by the truck envelope (fig3; midspan moment–curvature in fig2).
- Visible crack zone at rest at η = 0.36: 12.25 m (node spacing 0.61 m). On synthetic data this is a model output only (the masks are drawn from the same model); it becomes the spec's held-out check once field images exist.
- Mesh (fig6): the intact service deflection (stage-B dead load + truck) is 11.54 mm for N = 10…160 with a relative spread of 1e-9 (round-off): Hermite elements are exact at the nodes for uniform EI. Damaged (η = 0.36), relative difference to N = 160: 8.0e-4 (N = 10), 1.3e-7 (20), 5.9e-5 (40), 1.2e-5 (80), i.e. 0.0012 % from N = 80 to 160; not monotone, all below 0.1 %. The model uses N = 40.
- Secant loop: 1 pass at η = 1 (EC2) and at η = 0.36 (EC2 and Branson) — on the determinate girder the moments come from statics (P5 confirmed); the grillage (Plan 2) is where the loop matters.
- Branson interpolation (same EC2 section and Mcr) vs EC2 at η = 0.36: service deflection 73.3 vs 81.5 mm (Branson 10 % lower), truck-only 65.6 vs 72.1 mm, lowest EI ratio 0.13 vs 0.12.
- Monte Carlo stability (case2 at GSD 0.02), η p05 / p50 / p95: 0.338 / 0.360 / 0.385 (1,000 samples), 0.339 / 0.361 / 0.384 (10,000), 0.338 / 0.361 / 0.384 (2,000, seed + 1) — stable to about 0.001.
- Sensitivity (fig8; as-built 30-strand girder, case2 at GSD 0.02), swing of η over the tested range: girder spacing 1.8–2.6 m 0.082; strand centroid 80–120 mm 0.049; slab thickness 180–220 mm 0.029; bond factor ξ 0.5–0.7 0.027; strand cover 35–55 mm 0.014; girder fck 35–45 MPa 0.0009. The spacing swing is larger than the measurement band at GSD 0.02 (0.046 wide): at a fine GSD the stand-in geometry, not the camera, limits the estimate.
- Timing: linear solve 1.5 ms, nonlinear analysis 34 ms, inverse 2.7 ms, Monte Carlo (2,000) 5.2 s → no surrogate needed.
- ClaimBridge's own crack-width test suite passed (run locally on 2026-10-12).

**Limitations (for the report):**
- The synthetic masks are straight, constant-width bars: they test the pipeline link and its gates, not its measurement accuracy. The same model makes the widths and inverts them, so recovering η true shows consistency, not truth.
- Stand-in girder: AASHTO Type IV with a composite slab, strand position, concrete strengths and deck layout all stand-ins; the 1-D results use interior-girder rules (tributary width S, S/5.5), not the exterior girder G1.
- Not modelled: differential deck shrinkage and the timing of losses (all placed before casting); each is about as large as the staging shift (sizes in DECISIONS 2026-10-12, G09 review).
- The truck is placed for maximum midspan moment, which slightly underestimates the maximum live-load deflection (size in DECISIONS 2026-10-12); L/800 is compared with the truck deflection using S/5.5 and × 1.30, where AASHTO deflection checks often share all lanes equally.
- The pipeline's ± u is read as a Gaussian k = 2 band (σ = u/2), truncated at zero width; the η band covers measurement noise only.
- A crack seen at rest may be a residual opening from a past overload; η would then be underestimated (spec §5.6).

**Changed & why:** deviations from Plan 1, each in DECISIONS 2026-10-12: staged construction and HS20-44 × 1.30 (student decisions); pre-cracked zone from the moving-truck envelope; a third scenario (η_dec − 0.35) so that one case is measurable at GSD 0.10; the parser fails closed and an unresolved crack blocks a governing width; a refusal is reported as a lower bound on η (nominal and 95 %); η_ref is the first detectable case; the tornado holds the as-built 30 strands (re-designed variant kept as `tornado_redesigned`); deflection bands are stored in reversed order; figures rebuilt after review (fig1, fig2, fig6, fig7, fig8); the live run of `tests/test_claimbridge.py` now has 45 tests (plan: 5) after hardening, 3 of which call the private pipeline; the reading guide maps `f_eq` to the slides' {f0} with the same sign (the plan had f_eq = −F0).

**Numbers:** `results/stage1.json`; figures `results/figures/fig1…fig8.png`; `python -m pytest -q` → `256 passed, 3 skipped in 6.34s` (the 3 skips are the tests that call the private pipeline); with a pipeline config present, `python -m pytest tests/test_claimbridge.py -q` → `45 passed in 1.71s`.

**Reading guide:** for this block read uq.py → claimbridge.py → scripts/run_all.py; the course parts to defend in depth are:

| Code | Course | Where |
|---|---|---|
| `elements.beam_k`, `hermite` | beam element stiffness from Hermite shape functions | Logan (5th ed.) §4.1: shape functions p. 155 (Eq. 4.1.7), k p. 157 (Eq. 4.1.14) |
| `elements.udl_equivalent`, `point_equivalent` | work-equivalent nodal loads; `f_eq` is the slides' {f0} with the same sign ({f} = [k]{d} − {f0}); fixed-end reactions are −{f0} | Ch3 slides 3-8 to 3-10; Logan (5th ed.) §4.4 pp. 177–179 (Eqs. 4.4.3–4.4.11), Table D-1 p. 753 |
| `solve.solve` | assembly + partition K_ff D_f = F_f, reactions R = K D − F at the supports | Ch3 slides 3-6 to 3-7; Logan (5th ed.) §4.2–4.3 pp. 161–164 (Eqs. 4.2.3–4.2.5, 4.3.4, 4.3.7) |
| `nonlinear.solve_secant` | direct iteration with secant stiffness, D(i+1) = K_s(D(i))⁻¹ F — the same update as the slides' secant-stiffness algorithm | Ch11 slides, "Nonlinear Analysis Algorithms", 3) Secant Stiffness Algorithm |
| `rootfind.illinois` | secant method kept inside a bracket | Ch11 (secant idea) |

Explain-back questions: (1) Why does a 0.15 mm crack at rest mean prestress was lost, when the same crack in an RC beam is harmless? (2) Why can the field GSD (0.38 mm/px) not support this indicator, and what does a refusal at GSD 0.02 mm/px still tell you? (3) Why does the Monte Carlo band widen with GSD, and at GSD 0.02 mm/px what would narrow the η estimate more — a finer camera or the as-built drawings (compare fig5 with fig8)?

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 4 — B4–B5 Stage-2 grillage of S01 (2026-10-15)

**Planned:** put the Stage-1 girders into a grillage of the whole S01 deck, so that the FE model sits inside the secant loop and the η inverse (spec P5, §3 Stage 2); η from a crack on the exterior girder G1 and on the interior girder G3, with a Monte Carlo band; load shed to the neighbours and the deflection of every girder (Plan 2).

**Found:**
- Grid element: reproduces Logan (5th ed.) Example 5.5 — reduced stiffness, d1y = −2.83 in, φ1x = 0.0295 rad, φ1z = −0.0169 rad and every printed local force — to Logan's three figures, and the sympy derivation exactly; the same example through assembly and partition gives the same numbers.
- Linear checks: identical girders under equal load each carry qL²/8 and deflect 5qL⁴/384EI; a lone girder line is the Stage-1 girder (one pass); a torsion-free deck follows Courbon's rule with rigid cross members and the lever rule with flexible ones, and the real deck lies between them. PyNite 2.0.2 agrees with the grid solve (linear intact grillage, stage-B dead load + truck at G1's critical placement): max |Δw| = 2.66e-11 mm against a largest deflection of 11.89 mm; midspan moments to a relative 1.99e-12.
- Model as run: diaphragm = web, lane rule = AASHTO 12 ft design lanes, stage B in 10 truck increments, cracking flags one event at a time (DECISIONS 2026-10-15).
- Design: strands G1–G5 = 30, 29, 27, 29, 30 (Stage 1: 30) after 2 design rounds; soffit stress under service −0.435 / −0.439 / −0.326 MPa (G1 / G2 / G3, all ≤ 0). Truck midspan moment at each girder's critical placement 1652.8 / 1415.1 / 1189.2 kN·m (G1 / G2 / G3); two trucks side by side govern G1 (wheel lines at z = 0, 1.83, 3.66, 5.49 m: the outer wheel line over G1).
- Blind zone: no crack at rest above η_dec = 0.536 (G1) and 0.615 (G3); Stage-1 interior girder 0.557.
- η per measured case (Monte Carlo 5–50–95 %, 2,000 samples on the tabulated w(η); the exact inverse at the same widths differs by at most 1.4e-4; G1 case3 at GSD 0.10 drops 9 of 2,000 samples whose width exceeds what G1 can show at η = 0.02, every other case drops none). The synthetic cracks were generated with the Stage-1 single-girder (1-D) model, so the second column is the η used to generate the crack, not a truth for the grillage:

| Case | η used to generate the crack (1-D model) | Stage 1 | Grillage, G1 cracked | Grillage, G3 cracked |
|---|---|---|---|---|
| case2_gsd0p02 | 0.36 | 0.338–0.361–0.385 | 0.263–0.287–0.314 | 0.330–0.353–0.378 |
| case2_gsd0p05 | 0.36 | 0.312–0.367–0.433 | 0.236–0.294–0.374 | 0.304–0.359–0.434 |
| case3_gsd0p02 | 0.21 | 0.196–0.215–0.234 | 0.120–0.138–0.157 | 0.192–0.209–0.227 |
| case3_gsd0p05 | 0.21 | 0.169–0.214–0.262 | 0.093–0.138–0.185 | 0.166–0.209–0.255 |
| case3_gsd0p10 | 0.21 | 0.124–0.214–0.312 | 0.050–0.137–0.236 | 0.123–0.208–0.304 |

- The same measured crack implies a lower η in the grillage on G1: 0.287 against 0.361 on the Stage-1 girder (case2 at GSD 0.02), and 0.138 against 0.215 (case3 at GSD 0.02). Load shedding at rest reduces G1's moment: at η = 0.287 its stage-B dead-load midspan moment is 60.8 kN·m against 325.5 kN·m intact, so the girder needs less remaining prestress to open the same width (fig9: the G1 curve lies below the Stage-1 curve everywhere; the kinks in the grillage curves are discrete cracking events, one element flagged at a time). Stage A (prestress, girder, wet slab) comes from statics on the bare girder and stays on G1 whatever η: 1612.2 kN·m at midspan. Total midspan moment of G1 at rest (stage A + stage-B dead) 1937.8 kN·m intact → 1673.1 kN·m at η = 0.287, against 2002.2 kN·m on the Stage-1 girder at any η; in service (+ truck) 3590.5 → 2173.5 kN·m. G3 sheds at rest too (28.4 against 251.8 kN·m at η = 0.353) but stays close to Stage 1 (0.353 against 0.361); it also differs from the Stage-1 girder in its strands (27 against 30), so this near agreement is not a check of either model. At GSD 0.02 the structural model moves η (1-D → grillage, G1: 0.074) more than the measurement does (5–95 % band 0.051 wide), and the diaphragm and lane choices inside the grillage move it by up to 0.043 (sensitivity below).
- Refusals: nominal gate bounds η ≥ 0.384 / 0.295 / 0.195 at GSD 0.02 / 0.05 / 0.10 (G1) and 0.441 / 0.360 / 0.264 (G3); none at GSD 0.38.
- Indicator, G1 cracked at η = 0.287 (reference case case2 at GSD 0.02): stage-B service midspan moments G1–G5 intact 1978.3 / 1734.0 / 1328.3 / 1032.5 / 582.2 → damaged 561.2 / 2913.4 / 2012.8 / 1034.4 / 133.5 kN·m; G1 loses 71.6 % of its moment, G2 takes 59.6 % and G3 34.6 % of it, while G5's moment falls by 22.7 % of it (fig11). Truck-only deflection of G1 9.52 → 21.9 mm (band 21.1–21.9–22.2 mm), limit L/800 = 30.6 mm (fig12): inside the limit. Comparing this with the Stage-1 girder's 72.1 mm (Entry 3) is not like-for-like (different η, different girder); at the same η = 0.36 (`at_stage1_eta`) the truck-only deflection is 20.4 mm for G1 (intact 9.52) and 9.15 mm for G3 (intact 7.01) against 72.1 mm on the Stage-1 girder. The Stage-1 girder cannot shed: its moment comes from statics with a fixed distribution factor (S/5.5), so it keeps its whole truck share however soft it becomes; in the grillage the softened girder passes the truck to its neighbours (G1 at η = 0.36 loses 66.0 %, G2 cracks in 1 element). The verdict (inside L/800) holds for every diaphragm and lane variant at η = 0.287: G1 truck deflection 21.9 (web) / 23.9 (none) / 20.5 (composite) / 23.9 mm (clearance; 23.6 mm with its own strands). Cracked elements per girder 12 / 2 / 0 / 0 / 0 (G2 cracks under the shed load); secant updates in the last load step dead / service 6 / 16 (32 in the last service step with its re-converged cracking events).
- Indicator, G3 cracked at η = 0.353: moments 1088.4 / 1366.5 / 1441.0 / 1471.6 / 1287.8 → 1315.1 / 1779.7 / 161.3 / 1884.8 / 1514.4 kN·m; G3 loses 88.8 %, G2 / G4 take 28.7 / 28.7 % (G1 / G5 15.7 / 15.7 %); truck deflection of G3 7.01 → 9.16 mm (band 9.11–9.16–9.20 mm); cracked elements 0 / 0 / 12 / 0 / 0; updates 6 / 6.
- G1 model sensitivity at η = 0.287 (strands held at 30, 29, 27, 29, 30; reference case case2 at GSD 0.02):

| Diaphragm, lane rule | G1 loses (%) | w_pred (mm) | G2 cracked elements | η p50, reference case |
|---|---|---|---|---|
| web, lanes (as run) | 71.6 | 0.1604 | 2 | 0.287 |
| none, lanes | 61.6 | 0.1948 | 1 | 0.313 |
| composite, lanes | 76.0 | 0.1499 | 1 | 0.279 |
| web, clearance | 70.9 | 0.1604 | 2 | 0.287 |

  Own strand design per variant: "composite" keeps 30, 29, 27, 29, 30; "none" designs 30, 29, 28, 29, 30 (G3 + 1; G1 results unchanged: 61.6 %, 0.1948 mm, 1 element, η p50 0.313); "clearance" designs 31, 29, 27, 29, 31, and with its own strands G1 loses 69.5 %, w_pred = 0.1392 mm, G2 cracks in 2 elements and η p50 = 0.270 (η_dec 0.514). The diaphragm stand-in (S27) moves η more than the lane rule does.
- Secant loop on the grillage: more than one pass in every damaged analysis (at least 5 per load step; P5 shown on the indeterminate model), monotone decrease of the moment change (fig13); every analysis converged — all table points (a NaN would have stopped `build_table`), both indicators and all six sensitivity runs.
- Torsion sensitivity (GJ of cracked girder elements × 0.5): G1 w_pred 0.1604 → 0.1603 mm (−0.07 %), G1 stage-B service midspan moment 561.2 → 568.9 kN·m (+1.4 %); G3: no change at four significant figures (w_pred 0.1604 mm, moment 161.3 kN·m) — the interior girder is twisted little by the deck.
- Mesh (fig14), relative to N = 80: w_pred at N = 20 within 1.5e-5 (G1) and 4.6e-5 (G3); damaged service deflection at N = 20 within 1.2 % (G1: 26.67 against 26.36 mm; N = 10 0.57 %, N = 40 0.49 %, not monotone) and 4.9e-5 (G3); the model uses N = 20. The inverse quantity (w_pred) is mesh-converged; G1's damaged service deflection carries about 1 % mesh uncertainty (more than the plan's 0.6 % probe).
- Table: largest interpolation error at the knot midpoints 3.3e-3 (G1), 2.4e-3 (G3), relative to the largest tabulated width (so a smaller fraction of the widths actually measured); the meaningful check is at the η level: the Monte Carlo p05/p50/p95 on the table and the exact inverse at the same widths differ by at most 1.4e-4.
- Timing: grillage build 4.5 s, linear solve 11.4 ms, nonlinear service analysis 2.37 s, table 10.0 / 10.7 s (G1 / G3), all cases with Monte Carlo 32.0 / 32.8 s (G1 / G3), indicators at the Stage-1 η 5.8 s, G1 model sensitivity 124 s, mesh study 172 s; whole run 7.2 min.

**Limitations (for the report):**
- Every structural value is a stand-in (girder, strands, deck layout, overhang, diaphragms); the crack widths are synthetic and were generated with the Stage-1 single-girder model, so the grillage bands cannot be scored against a truth: they show how far the structural model moves η.
- The diaphragms are their webs alone, rigidly connected to the girder nodes: a lower bound as a beam, an upper-bound assumption at the connection, so the model does not bound the load sharing; the three-case sensitivity above shows the spread (G1 loses 61.6–76.0 %). Bearings leave the girder ends free to twist.
- The initial pre-cracked zone comes from a linear intact envelope (the cracked girder's own shedding is ignored there, so it errs toward more pre-cracking); the truck stays at its Stage-1 longitudinal placement and at the critical transverse placement of the intact deck.
- The overhang's wet-slab torque on the bare exterior girder (stage A) is not modelled; stage A is vertical statics, as in Stage 1. Stage A from statics also means that prestress lost after composite action is not restrained by the deck (the neighbours take back none of the lost prestress effect); modelling that would push the G1 η estimate lower still, the same direction as the shedding.
- No tornado on the grillage (DECISIONS 2026-10-15); the Stage-1 tornado stands.

**Changed & why:** deviations from Plan 2, each in DECISIONS 2026-10-15 (Plan 2 group B fixes and the close-out): lane rule = AASHTO design lanes, incremental stage-B loading, cracking flags one event at a time, diaphragm three-case option (all before this block); in this block the η that generated the synthetic crack is labelled `eta_gen_1d` ("η used to generate the crack (1-D model)") instead of "η true"; a G1 model-sensitivity block (diaphragm none / web / composite, lane rule clearance, own strand design per variant) was added to `run_stage2.py` and `stage2.json`; fig12 shows the truck-only deflection against the L/800 line instead of the dead + truck deflection; `run_stage2.py --figures-only` redraws fig9–fig14 from `stage2.json`; after the final review, `at_stage1_eta` (G1 and G3 at the Stage-1 η = 0.36, a like-for-like comparison), total midspan moments (stage A + stage B) and truck deflections per sensitivity row were added to `stage2.json`, and the fig11 title says it shows stage B only; the whole run takes about 7 minutes (plan: 6; the sensitivity block adds about 2).

**Numbers:** `results/stage2.json`; figures `results/figures/fig9…fig14`; `python -m pytest -q` → `401 passed in 59.64s` with the pipeline configuration (`398 passed, 3 skipped in 58.44s` without it).

**Reading guide:** read elements.py (grid part) → gridmodel.py → gridsolve.py → nonlinear.secant_iterate (the course part), then deck.py → grillage.py → wtable.py → scripts/run_stage2.py. The parts to defend in depth:

| Code | Course | Where |
|---|---|---|
| `elements.grid_k_local`, `grid_T`, `grid_k_global` | grid element: beam bending + torsion bar, transformation k = Tᵀ k̂ T | Logan (5th ed.) §5.4: Eq. 5.4.16 p. 241, Eqs. 5.4.17–5.4.18 p. 243, Eq. 5.4.19 p. 244; Example 5.5 pp. 244–250 |
| `gridsolve.stiffness`, `solve` | assembly member by member, partition K_ff D_f = F_f, reactions, local end forces f̂ = k̂ T d − f̂_eq | Ch3 slides 3-6 to 3-10; Logan §5.4 |
| `nonlinear.secant_iterate` | direct iteration with secant stiffness on an indeterminate model, D(i+1) = K_s(D(i))⁻¹ F | Ch11 slides, "Nonlinear Analysis Algorithms", 3) Secant Stiffness Algorithm |
| `deck.lever_weights`, `area_to_girders` | statically equivalent load transfer to the girder lines | statics (lever rule) |

Explain-back questions: (1) Why does the secant loop need more than one pass on the grillage when it needed exactly one on the single girder? (2) The same crack width gives a different η on G1, on G3 and on the Stage-1 girder (fig9, fig10): which model differences cause that? (3) When G1 cracks, G2 takes part of its moment: what does that do to G2's zero-tension design, and what would you check on the real bridge?

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 5 — B6 submission packaging and report QA (2026-10-18)

**Planned:** convert the technically complete repository into a reviewable submission candidate: close factual placeholders, expose the important results outside the large JSON files, capture a reproducibility record, and produce matched DOCX/PDF files.

**Found:** the implementation was ahead of the original timetable, but the deliverable was not yet closed. The report still named an old results commit, contained three student placeholders, described the test count ambiguously, and existed only as a draft DOCX plus an internal preview PDF. The progress log recorded extensive technical work but no completed student review.

**Changed & why:** the factual placeholders were replaced without inventing student experience; the ClaimBridge boundary now states that the course repository begins at the mask-and-GSD interface; the discussion labels the result as a screening indicator; the report points to a separate student review record; `results/key_results.csv` makes the Stage-1/G1/G3 inverse comparable without opening JSON; `results/run_record.json` records the result provenance; and the submission candidate is rendered as matched DOCX and PDF files. A workflow figure and a submission README make the project understandable without reading the implementation plan.

**Numbers:** `python -m pytest -q --basetemp <workspace temp>` → 401 passed in 108.18 s with the configured private pipeline. The result files remain the committed Stage-1 and Stage-2 outputs from 2026-10-15; this packaging block does not recompute or change them.

**Reading guide:** start with `SUBMISSION_README.md`, then `results/key_results.csv`, then the report. Use `results/run_record.json` to identify the code and data snapshot. Complete `report/STUDENT_REVIEW.md` before changing the report status from submission candidate to final.

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 6 — ClaimBridge supervisor background note (2026-10-18)

**Planned:** make the submission understandable to an examiner who did not participate in the ClaimBridge thesis work, without exposing the private pipeline or overstating the present validation.

**Changed & why:** added a standalone ClaimBridge background and FEM interface note in Markdown, Word and PDF. It defines the private/public boundary; documents the mask, surface-scale, resolution-gate and parsed-output contract; explains GSD and the three-pixel rule; follows one synthetic case from an eight-pixel mask to `0.16 ± 0.04 mm` and then to the Stage-1/G1/G3 inverse results; separates measurement uncertainty from model uncertainty; and adds three peer-reviewed image-based crack-quantification references for context. The note states that these publications contextualise the measurement problem but do not validate ClaimBridge.

**Evidence:** `docs/CLAIMBRIDGE_BACKGROUND.md`, `report/ClaimBridge_Background_and_FEM_Interface.docx`, `report/ClaimBridge_Background_and_FEM_Interface.pdf`, `src/crackfem/claimbridge.py`, `src/crackfem/uq.py`, `data/claimbridge/case2_gsd0p02_spike.json`, and `results/key_results.csv`.

**Reading guide:** read the Main point, Sections 2–6, and the Recommended interpretation. Be able to state which part is AI measurement, which part is FEM inference, why the GSD 0.38 mm/px case is refused, and why one width produces different eta values in the single-girder and grillage models.

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 7 — Project lineage and FEM deliverable boundary (2026-10-18)

**Planned:** make the narrative explicit for a supervisor who knows the FEM assignment but did not work on the thesis: identify ClaimBridge as prior background, state the hand-over variable, and list the new work delivered to the FEM course.

**Changed & why:** revised the background note, report introduction and submission README around one continuous chain: `ClaimBridge background → w ± u → inverse FEM → η → nonlinear girder/grillage response → engineering outputs`. The documents now state that ClaimBridge supplies the measurement observation and is not being resubmitted as the FEM contribution. A four-row deliverable matrix identifies the verified Stage-1 beam model, inverse structural interpretation, Stage-2 nonlinear grillage and reproducible evidence package as the course outputs.

**Evidence:** `docs/CLAIMBRIDGE_BACKGROUND.md`, report Section 1.2, `SUBMISSION_README.md`, the Stage-1 and Stage-2 result files, 401-test record and submission artifacts.

**Reading guide:** explain the project in three sentences: what existed before the course, what engineering gap remained, and what new FEM outputs close that gap.

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 8 — B6 closing the proposal's outcomes (2026-10-24)

**Planned:** close the four Priority-1 items of the 2026-10-18 compliance audit that do not need the student: a compliance table in the report, the width-response curves of Expected Outcome 2, a number for Expected Outcome 3, and a visible hand calculation. Commit the uncommitted 2026-10-18 work first and journal it.

**Found:** the 2026-10-18 work had been left uncommitted on branch `report` with no DECISIONS lines (committed as f3b6ce5, journal e0ef04b). The candidate report had been built before the audit, so none of the items was in it. The proposal's one-element damage model gives an answer set by the mesh: at η = 0.36 its truck deflection is 27.1 / 18.7 / 14.2 mm at N = 10 / 20 / 40 against 72.1 mm for the distributed model (28 / 15 / 7.5 % of the damage-induced deflection), and at N = 10 it would have passed L/800. On the Stage-1 girder the truck deflection already reaches 38.0 mm at η_dec, where no crack is open at rest.

**Changed & why:** `Girder1D` got an optional `damaged` mask (elements outside keep their uncracked EI; default unchanged) and `damaged_zone(length)`; `scripts/run_proposal_checks.py` writes `results/proposal_checks.json` and fig15–fig17 without touching stage1.json or stage2.json; `figures.fig_width_response` and `fig_damaged_length` are new. The report has a new Section 8 (Block 6: 8.1 curves, 8.2 damaged zone, 8.3 hand calculation, 8.4 compliance table) and Appendix B; old Sections 8–13 are now 9–14. `scripts/build_submission_report.py` renders the new markdown sections into the Word template (headings, bullets, numbered items, tables, figures) so the text is written once.

**Numbers:** `results/proposal_checks.json` (run 164 s); `python -m pytest -q` → 416 passed with the pipeline configuration. Hand calculation: midspan 2.0338 mm = 5qL⁴/384EI, reactions 52.947 kN each. Width response at w = 0.16 ± 0.04 mm: Stage-1 deflection 69.1 – 72.1 – 74.4 mm, α 0.128 – 0.132 – 0.138, midspan EI ratio 0.119 – 0.122 – 0.126; all ratios 0.100 – 0.113 above w = 0.4 mm.

**Reading guide:** girder1d.py (`damaged_zone`, the `damaged` branch in `stiffness_fn`) → tests/test_damaged_zone.py → scripts/run_proposal_checks.py (`hand_calc` first, then `stage1_rows`, `damaged_length_study`) → report Section 8 and Appendix B. Redo Appendix B by hand on paper (one element, then node 2) before reading the code's numbers.

Explain-back questions: (1) Why does the one-element result fall toward the intact deflection as N grows, while the fixed 4.90 m zone does not change? (2) Why can the Stage-1 girder exceed L/800 at η_dec although no crack is open at rest, and why does the grillage girder not? (3) Why is the equivalent uniform α given for the single girder only?

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 9 — Closing block before the hand-in (2026-11-01)

**Planned:** close the last audit items that do not need the student, check the report against the hand-in requirements (more than 20 pages, progress explained clearly, a clean layout), and prepare a public repository with the history squashed into one commit. The student review was moved to the end at the student's request.

**Found:** audit item 6 (the proposal's w ± u end points) was still only a sentence. The PDF was 23 pages without a cover or contents page; the design-review and other draft tables had equal-width columns (the "#" column as wide as the text columns), rows split across pages, and two block labels had lost their bold when the builder rewrote the paragraph. The report META line cited commit 816b300, which will not exist in the squashed public history, and the Word file carried the e-mail address of the Word account as "last modified by".

**Changed & why:** `end_point_rows` in `scripts/run_proposal_checks.py` solves η again at w − u and w + u for the five measured Stage-1 cases, and `main` checks that η at w matches stage1.json and that the end points enclose the Monte Carlo band; the table is in report Section 8.1 and the compliance row now says "Delivered". `build_submission_report.py` adds a cover page, a contents page (a TOC field that Word fills when the PDF is exported), the "Progress at a glance" table of Section 1.3 (rendered from report_draft.md), content-based column widths, rows kept whole with repeated headers, bold block labels and appendices on a new page; the META line gives the repository URL instead of a commit hash. README gained a "History" note on the squash.

**Numbers:** `results/proposal_checks.json` rerun (130 s), unchanged apart from the new `end_points`; figures 15–17 bit-identical. End points at w = 0.16 ± 0.04 mm (GSD 0.02): η = 0.334 – 0.390 against the Monte Carlo 0.338 – 0.385; over the five cases the end points are 19–22 % wider than the band. `python -m pytest -q` → 420 passed with the pipeline configuration. Report: 26 pages. xray: only the {f0} vector notation flagged as template variables and three em dashes inside cited titles (`report/xray/xray_2026-11-01_closing.json`).

**Reading guide:** `end_point_rows` and the two checks in `main` (scripts/run_proposal_checks.py) → the end-point table in report Section 8.1 → the Progress at a glance table in Section 1.3.

Explain-back questions: (1) Why must η at w + u be below the 5 % point of the Monte Carlo band? (2) Why is the end-point interval about 2/1.645 times as wide as the 5–95 % band? (3) Why can a squashed public history not keep the commit hashes quoted in the log?

**Reviewed by student:** declared complete by the student, 1 November 2026.

## Entry 10 — Browser demo for the instructor (2026-11-01)

**Planned:** a link the instructor can open and try without installing anything: set a crack width, see η, the stiffness along the span and the truck deflection.

**Found:** the saved results already hold the width–response rows of every model (figure 15), but the grillage was solved at only 14 values of η per girder, too coarse to move a slider through. One grillage point near decompression (G3, η ≈ 0.595) did not converge with the default secant update.

**Changed & why:** `scripts/build_demo_data.py` re-solves the single girder at 81 values of η and the grillage at 31 per cracked girder (G1, G3) and writes `docs/demo_data.js`; a grillage point that does not converge is retried once with relax 0.5 and otherwise dropped and listed (none was dropped). `docs/index.html` is one page with no library: a crack-width slider and a GSD choice, three result cards (η with the w ± u band, truck deflection against L/800), the stiffness ratio of every element along the span, the width–deflection curves with the refused strip and the blind zone, and the stage-B moment of each girder. It interpolates between solved points and runs no FEM. GitHub Pages serves it from `main`/`docs`. `tests/test_demo_data.py` checks the page's interpolation against the saved results.

**Numbers:** at w = 0.16 mm and GSD 0.02 the page gives η 0.361 (band 0.334–0.390), 0.287 (G1) and 0.353 (G3), truck deflection 72.0, 21.8 and 9.2 mm; a refused 0.03 mm reading gives η ≥ 0.439, as `stage1.json` (case1). Data export 322 s. `python -m pytest -q` → 424 passed (420 + 4 new). Checked in the browser at desktop and phone width, light and dark. xray on the page text: clean.

**Reading guide:** `atWidth` and `evaluate` in `docs/index.html` (the interpolation and the band) → `eta_at` in `tests/test_demo_data.py` (the same interpolation in Python) → `grillage_models` in `scripts/build_demo_data.py`.

Explain-back questions: (1) Why can the page show a deflection band without running the FEM? (2) Why is the band of a refused reading one-sided? (3) Why does the curve jump at zero width?

**Reviewed by student:** declared complete by the student, 1 November 2026.
