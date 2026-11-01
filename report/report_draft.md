% TITLE: From a Measured Crack Width to Remaining Prestress, Stiffness and Load Sharing
% SUBTITLE: An FEM–AI structural indicator for a prestressed girder bridge (CPK72, span S01)
% META: 01203541 Finite Element Analysis of Structures · Kasetsart University · 1st semester 2026
% META: Course project (solo, 15 %, CLO9 "AI + FEM") · Progress report
% META: Palid Archavachalee (6914500959)
% META: 1 November 2026 · code, tests and results: https://github.com/Palarch/crackfem

# Summary

ClaimBridge existed before this course project as part of the author's thesis work on image-based bridge inspection. It supplies a measured crack width and a resolution band; it does not calculate the structural consequences. The new FEM course work begins at that boundary and maps the ClaimBridge measurement to a model-conditional estimate of remaining prestress (η, a fraction of the design value), stiffness, truck deflection and load sharing.

The finite element model is the engine of that answer. Each girder element gets its stiffness from its own moment through the Eurocode 2 cracked-section rules, so the stiffness depends on the solution; on the five-girder grillage of the deck the moments depend on every girder's stiffness, so the answer needs a nonlinear (secant) FE solve inside an inverse loop for η.

The report follows the work block by block, as it happened. The approved proposal had five problems that changed the design (Block 0). A one-girder model (Stage 1) was built and verified first, then the whole deck as a grillage (Stage 2). Main findings, all on a stand-in girder and on synthetic crack images:

- A girder can lose up to 44 % of its prestress without showing any crack at rest (η above 0.557 gives no crack), so a crack-based indicator is blind there.
- At the pipeline's field resolution (0.38 mm per pixel) no crack this girder can show is measurable; at 0.02 mm per pixel a 0.16 mm crack gives η = 0.338–0.385 (5–95 %).
- The rule-engine rating is 5 "Excellent" for every measured case, including a crack that the stand-in model maps to a 79 % prestress loss: the rating carries no information about prestress, which is the gap the indicator addresses.
- The structural model matters more than the camera at fine resolution: the same crack gives η = 0.361 on the single girder and 0.287 on the exterior girder of the grillage, and at the same η the single girder deflects 72.1 mm under the truck (2.4 times the L/800 limit of 30.6 mm) where the grillage girder deflects 20.4 mm (inside the limit), because the grillage girder sheds about two thirds of its load to its neighbours.
- Putting all the damage in one element, as the proposal planned, makes the answer depend on the mesh: at η = 0.36 one element gives 27.1 mm of truck deflection at N = 10 (inside L/800) and 14.2 mm at N = 40, where the model with a cracked zone set by the moments gives 72.1 mm at every mesh.

The code (Python, 425 automated tests) was written by Claude (Anthropic's AI). The student's final code and oral review is recorded separately before submission; Section 1.4 states the arrangement.

# 1 Introduction

## 1.1 The problem

Thai highway bridges are inspected visually and each defect gets a condition rating. A crack is rated by its width. In a reinforced concrete beam a hairline flexural crack is expected and may be acceptable. In a prestressed girder designed to remain in compression under permanent load, an open flexural crack at rest is evidence that the local precompression may have been exhausted. Under the adopted model, and provided the opening is not residual damage from a past overload, its width carries information about the prestress that remains; a categorical rating does not use that information.

## 1.2 From ClaimBridge background to the FEM deliverables

The project follows one continuous chain: ClaimBridge thesis background → crack width w ± u → inverse FEM → remaining prestress η → nonlinear girder and grillage response → engineering outputs. The hand-over variable is the measured width and its resolution band. ClaimBridge is the upstream measurement source; the structural-analysis layer is the work submitted to the FEM course.

- **AI part (measurement).** ClaimBridge, the author's thesis inspection pipeline (private), receives a crack mask and surface scale, measures the crack width, and passes the result to a rule engine. The upstream image-to-mask AI remains outside this course repository and is not evaluated here; the interface tested by this project begins with the mask and GSD. The pipeline reports the width as w ± u with u = 2·GSD and refuses cracks narrower than 3 pixels.
- **FEM part (mechanics).** The finite element model of the girder or the whole deck maps that width to a model-conditional estimate of the remaining prestress η (an inverse problem), and η to stiffness, deflection and load sharing (a forward nonlinear analysis).

The FEM course deliverables are a verified Stage-1 Euler–Bernoulli beam model, an inverse solver for η with a Monte Carlo measurement band, a nonlinear five-girder grillage that captures load redistribution, and a reproducible evidence package containing tests, result files, figures, sources and decision records. ClaimBridge is not resubmitted as the FEM contribution; it supplies the observation that makes the FEM problem possible.

The result is an **advisory structural indicator** placed next to the pipeline's rating. It never changes the rating.

![Figure 0. Boundary and data flow from the crack image to the engineering outputs. The private upstream image-to-mask AI is outside the course repository.](../results/figures/fig0_workflow.png)

## 1.3 How to read this report

The course asked for a report that tells the progress as it happened. Sections 3 to 8 follow the work blocks in order. Each block says what was planned, what was found, what changed and why, and the numbers it produced. Every number comes from the result files in the repository (`results/stage1.json`, `results/stage2.json`) or from test output, and every non-obvious decision has a dated line in `docs/DECISIONS.md`. Section 8.4 maps every promise of the proposal to where it was delivered. Section 9 shows where the course's FE methods sit in the code, Section 10 collects the verification, and Sections 11 to 13 discuss the results, the limitations and the next steps. The separate ClaimBridge background note documents the inherited measurement side for a reader who was not involved in the thesis.

**Progress at a glance.** One row per work block. Each row is a section of this report and an entry in `docs/PROGRESS_LOG.md`; the test count is the automated test suite at the end of the block.

| Block | Dates | Planned | Delivered | Tests | Section |
|---|---|---|---|---|---|
| 0 Design review | 6–7 Oct | turn the approved proposal into a design | five problems found (P1–P5); the unknown became η; a staged plan (one girder, then the deck) | – | 3 |
| 1 Source pack | 7 Oct | a source for every input value | rows S01–S30 tagged primary, secondary or stand-in; three unsourced values replaced | – | 4 |
| 2 Stage 1 core | 7–12 Oct | one girder: section, crack width, stiffness, beam FE, inverse | checked against two published examples; the one-pass secant loop proved (P5) | 165 | 5 |
| 3 Pipeline link | 12–13 Oct | synthetic cracks through the real pipeline; Monte Carlo band | blind zone above η = 0.557; nothing measurable at GSD 0.38 mm/px; η = 0.338–0.385 at 0.02 | 256 | 6 |
| 4–5 Grillage | 15 Oct | the five-girder deck, a real secant loop, the inverse on G1 and G3 | the same crack gives η = 0.287 on G1; the L/800 verdict reverses | 401 | 7 |
| 6 Proposal outcomes | 18–24 Oct | close every promise of the proposal | width–response curves, one element against a damaged zone, hand calculation, compliance table | 416 | 8 |
| Closing | 1 Nov | the last audit items and the final build | w ± u end points beside the Monte Carlo band; report, repository, browser demo and review record | 425 | 8.1, 14 |

## 1.4 Authorship

The code was written by Claude (Anthropic's AI); this was the student's choice at the start of the design (DECISIONS 2026-10-07), because the final exam fell in the middle of the project. The student made the engineering decisions recorded as "student decision" in the log, including staged construction and the DOH truck factor. Each block in `docs/PROGRESS_LOG.md` ends with a reading guide, three explain-back questions and a `Reviewed by student` line. The student declared the review complete on 1 November 2026 (Section 14). The report text was drafted by Claude from the progress log and the decision journal; the decisions marked as the student's are the student's.

# 2 Starting point: the approved proposal

The proposal "FEM–AI Structural Indicator for Crack-Width-Informed Deflection: A Case Study on Bridge CPK72 (Route 304)" (15 September 2026) was approved without comments. It planned:

- a 1-D Euler–Bernoulli beam of span S01 (24.50 m, simply supported), modelled as an equivalent reinforced concrete rectangle;
- the crack width from ClaimBridge setting the effective moment of inertia of one cracked element through Branson's equation;
- an iteration between the effective stiffness and the applied moment until the deflection stops changing;
- a mesh convergence check with 10, 20 and 40 elements against the closed-form deflection;
- a ±2·GSD band on the result from the pipeline's measurement uncertainty, and a machine-learning surrogate as future work.

Section 8.4 maps each of these items and the proposal's five Expected Outcomes to the place where it was delivered.

# 3 Block 0: Design review (6–7 October)

**Planned:** turn the proposal into a design and an implementation plan.

**Found:** five problems, each of which would have made the result either wrong or not a real use of FEM.

| # | Proposal said | Problem found | Change |
|---|---|---|---|
| P1 | The crack width sets the effective stiffness through Branson's equation | Branson's equation has no crack-width term; it uses the moment ratio Ma/Mcr. The width → stiffness link would have been an invented curve | An explicit chain: crack width → strand stress increase (inverted EC2 §7.3.4 crack-width formula) → section state → stiffness |
| P2 | Damage in one element; mesh check with 10/20/40 elements against the closed form | The damaged length is L/N, so the result depends on the mesh (2.45 m against 0.61 m). And Hermite beam elements with work-equivalent loads give exact nodal deflections for any N, so matching the closed form proves nothing | Each element gets its stiffness from its own moment, so the cracked zone appears wherever M > Mcr. Convergence is measured on the damaged model |
| P3 | Equivalent RC rectangle 0.45 × 1.70 m | A 24.5 m highway span is typically a precast prestressed I-girder (the as-built drawings were not available) | Stand-in standard prestressed I-girder (AASHTO Type IV) with a composite slab, prestress sized by a design rule |
| P4 | The crack width comes straight from ClaimBridge | The pipeline has not yet measured a real crack on this bridge (first field survey: November 2026) | Synthetic crack masks of known width run through the real crack-width command; stated plainly throughout |
| P5 | Iterate effective stiffness ↔ moment until the deflection stops changing | S01 is simply supported, so the moments come from statics and never change: the loop ends after one pass | Model the whole deck as a grillage, which is statically indeterminate: the loop becomes a real nonlinear solve. The 1-D stage keeps a test that proves the one-pass behaviour |

**Changed and why:** the unknown became the remaining prestress η because, under the adopted flexural model, a crack visible at rest indicates that local precompression has been exhausted; η is therefore the model quantity inferred from the width. Two other approaches were rejected: a forward-only chain from width to moment to stiffness (no real nonlinearity), and continuous spans (no evidence that the spans are continuous). The work was split into Stage 1, one girder with the whole chain, complete and submittable on its own, and Stage 2, the grillage of the deck as the main deliverable. Crack width and stiffness both use Eurocode 2 (one code family, and EC2 has an explicit crack-width formula for prestressed members that can be inverted in closed form); ACI Branson is kept as a cross-check. The machine-learning surrogate was dropped, because one solve costs milliseconds (Section 6 gives the timing).

**Numbers:** a probe on 7 October confirmed that the pipeline's crack-width command measures a synthetic 6-pixel crack at GSD 0.05 mm/px as 0.30 ± 0.10 mm, and that the rule engine rates it 5 (Excellent).

# 4 Block 1: Source pack (7 October)

**Planned:** open every source behind the input values and lock them.

**Found:** every value got one of three tags: [primary] (the original was opened), [secondary] (an open document that quotes the original) or [stand-in] (a chosen value, not the real bridge). The paywalled standards (EN 1992-1-1, ACI 318, the AASHTO Standard Specifications) were read through open documents that quote them: JRC and European Concrete Platform worked examples, the Sétra bridge guide, TxDOT/TTI design examples, PCI Journal papers and NCHRP reports. The full table is `docs/sources.md` (rows S01–S30). The stand-in structure:

| Item | Value | Tag |
|---|---|---|
| Span S01 | 24.50 m, simply supported | [primary] DOH as-built drawing |
| Deck | 5 girders at 2.20 m, slab 0.20 m, wearing surface 0.05 m, 2 barriers | [stand-in] |
| Girder | AASHTO Type IV (1.37 m deep), composite slab | [primary] outline, [stand-in] choice |
| Strands | 12.7 mm, 1860 MPa, Ep = 195 GPa, jacking 0.75 fpu, lump-sum loss 310 MPa | [secondary] |
| Concrete | girder fck 40 MPa, slab 28 MPa, EC2 Table 3.1 for fctm and Ecm | [stand-in] / [secondary] |
| Crack width | EC2 §7.3.4, k1 = 1.6 (tendons), k2 = 0.5, k3 = 3.4, k4 = 0.425, kt = 0.4, bond factor ξ = 0.6 | [secondary]; ξ unverified |
| Stiffness | EC2 §7.4.3, β = 0.5 (sustained or repeated load) | [secondary] |
| Live load | HS20-44 × 1.30 (DOH standard), impact 50/(L + 125) ≤ 0.30 | [secondary] |
| Limit | live-load deflection L/800 = 30.6 mm | [secondary] |

- Three published values of the Type IV second moment of area disagree slightly; the exact integral of the outline in the code matches the NCDOT table to the last digit.
- The bond factor ξ = 0.6 for pretensioned strand could not be found in any open source; it is marked unverified and appears in the sensitivity study (Section 6).
- Three families of unit weights disagree (AASHTO, EN 1991-1-1, JRA). AASHTO was adopted so that dead and live load come from the same family.

**Changed and why:** three values without a source were replaced by sourced ones: barrier weight 7.0 → 4.76 kN/m (TxDOT T501 barrier), concrete unit weight 24.0 → 23.56 kN/m³ and wearing surface 22.0 → 21.99 kN/m³ (AASHTO 150 and 140 pcf).

# 5 Block 2: Stage 1 core, one girder (7–12 October)

**Planned:** build the section model, crack width, stiffness, beam FE, nonlinear loop and inverse, each with tests against independent evidence.

## 5.1 The physics chain

1. **Section.** An elastic fibre model of the composite girder: concrete linear in compression and, once cracked, zero in tension; strands bonded and elastic, carrying a locked-in strain from the prestress (η times the design value). For a given moment it returns the curvature, the strand stress increase Δσp and the cracking moment Mcr(η).
2. **Staged construction** (student decision, 12 October). Stage A (prestress, girder and wet slab) acts on the bare girder; stage B (wearing surface, barriers and truck) acts on the composite girder. The slab carries no stress from its own weight; it carries the stage-A strain plane as a locked-in strain.
3. **Crack width**, EC2 §7.3.4 written for a prestressed member, with Δσp in place of the reinforcement stress:
   - wk = sr,max · (εsm − εcm)
   - sr,max = k3·c + k1·k2·k4·φ / ρp,eff, with ρp,eff = ξ·Ap / Ac,eff
   - εsm − εcm = the larger of [Δσp − kt·(fct,eff / ρp,eff)·(1 + αe·ρp,eff)] / Ep and 0.6·Δσp / Ep
   The width is linear in Δσp on each branch, so the inverse w → Δσp is closed form. For every realistic width the 0.6 floor governs.
4. **Stiffness**, EC2 §7.4.3: curvature κ = ζ·κcracked + (1 − ζ)·κuncracked with ζ = 1 − β·(Mcr / M)² and β = 0.5. The FE model carries the stage-B loads, so each element's secant stiffness is measured from the end of stage A: EI = (M − MA) / (κ(M) − κ(MA)).

## 5.2 The FE part

- **Beam element** (course Ch. 3; Logan §4.1): Hermite shape functions, 2 DOF per node (deflection, rotation), 4 × 4 stiffness matrix; distributed and point loads as work-equivalent nodal loads {f0}.
- **Assembly and partition** (Ch. 3; Logan §4.2–4.3): global K assembled element by element, the free block solved K_ff·D_f = F_f, reactions R = K·D − F at the supports.
- **Secant iteration** (Ch. 11, "Secant stiffness algorithm"): solve with the current element stiffnesses, recompute each element's moment and secant EI, solve again, D(i+1) = Ks(D(i))⁻¹·F, until both the moments and the stiffnesses stop changing.
- **Inverse:** find η such that the predicted crack width at rest w_pred(η) equals the measured one. Illinois method (a secant method that always keeps the root bracketed), with bisection for the edge where the crack closes; a width below the cracking threshold returns the bound "η ≥ η_dec", never a point value.
- The 40-element girder is the default mesh.

## 5.3 What was found

- Beam nodal deflections are exact for N = 1, 2, 3, 5, 8 and 40 elements: the proposal's mesh study on the uncracked beam would have proved nothing (P2 confirmed).
- The secant loop stops after one pass on the simply supported girder, with or without damage (P5 confirmed). On an indeterminate test beam (a propped cantilever with softening stiffness) it needs at least two passes and moves moment from the fixed end to the span: there the loop is real. This is why Stage 2 exists.
- The model reproduces two published examples: the EC2 crack-width example of the European Concrete Platform (Example 7.3, wk = 0.184 mm) and a cracked prestressed section by Mast (PCI Journal 1998, Example 1: neutral-axis depth 17.26 in, strand stress 171.96 against 171.97 ksi printed). The staged section matches a two-stage transformed-section hand calculation within 1e-5.
- A section with its strands below the kern cracks at the top under prestress alone, because in cracked mode concrete carries compression only everywhere; the plan's test assumed otherwise and was corrected.

## 5.4 What the reviews caught

The two-stage reviews after each task group found real errors in the plan, not only in the code:

- a test compared the cracked and uncracked states at a moment where they cannot agree (fixed at 300 kN·m);
- a round-off bound of 1e-12 was below the noise floor (now 1e-9), and a moment reaction was summed into a vertical force check;
- the plan wrote the work-equivalent loads as f_eq = −{F0}; checking the Ch. 3 slides showed that the slides' {f0} is the equivalent load itself, so f_eq = +{f0} and the fixed-end reactions are −{f0};
- the size of the staging effect was first quoted as 0.15 in η; the measured shift is +0.06 to +0.10 (the first estimate mixed two inconsistent assumptions), and the student was told;
- the truck scan on a 1 mm grid could miss the moment peak, so the exact critical axle positions were added; the secant loop now also checks the stiffness change, because a moment-only check can stop with a stiffness that belongs to other moments where EI jumps at Mcr.

**Numbers:** 165 tests passed at the end of the block.

# 6 Block 3: Link to the thesis pipeline and Stage 1 results (12–13 October)

**Planned:** run synthetic cracks through the real pipeline, invert for η with a Monte Carlo band, check mesh convergence and the Branson cross-check.

## 6.1 Synthetic cracks

Three crack scenarios were generated with the Stage-1 model at η = η_dec − 0.10, − 0.20 and − 0.35 (η = 0.46, 0.36 and 0.21), drawn as straight, constant-width masks and measured by the pipeline at four GSDs (0.02, 0.05, 0.10 and 0.38 mm/px): 12 cases. Only files cross the boundary between the two code bases (a mask image and a surface file go in, a text report comes out); on 13 October the course repository was cleaned of every internal name of the private pipeline, and the command that runs it is configured per machine. The Monte Carlo draws 2,000 widths from a normal distribution with σ = u/2 (the pipeline's ±u read as a 95 % band), truncated at zero, and inverts each one.

## 6.2 Results

- **Design:** the smallest strand count with no tension at the soffit under service at η = 1 is 30 strands (P0 = 3212 kN).
- **Blind zone:** no crack is predicted at rest for η ≥ η_dec = 0.557. A girder can lose up to 44 % of its prestress without a visible sign. Below η_dec the predicted width grows as η falls, up to 0.89 mm at η = 0.02 (Figure 1).
- **Resolution:** the pipeline refuses cracks narrower than 3 px. At GSD 0.38 mm/px (the lower end of its field GSD) the gate is 1.14 mm, wider than any crack this girder can show (0.89 mm), so nothing is measurable and a refusal says nothing. At finer GSD a refusal still bounds η from below: η ≥ 0.439 / 0.368 / 0.272 at GSD 0.02 / 0.05 / 0.10 mm/px. In all, 7 of 12 cases were refused (Figure 2).

| Case | GSD (mm/px) | Width (px) | w ± u (mm) | η, 5 – 50 – 95 % | η used to generate | Rating |
|---|---|---|---|---|---|---|
| case2 | 0.02 | 8 | 0.16 ± 0.04 | 0.338 – 0.361 – 0.385 | 0.36 | 5 |
| case2 | 0.05 | 3 | 0.15 ± 0.10 | 0.312 – 0.367 – 0.433 | 0.36 | 5 |
| case3 | 0.02 | 20 | 0.40 ± 0.04 | 0.196 – 0.215 – 0.234 | 0.21 | 5 |
| case3 | 0.05 | 8 | 0.40 ± 0.10 | 0.169 – 0.214 – 0.262 | 0.21 | 5 |
| case3 | 0.10 | 4 | 0.40 ± 0.20 | 0.124 – 0.214 – 0.312 | 0.21 | 5 |

- The generating η lies inside every band, and the band widens with GSD because u = 2·GSD. The bands are measurement noise only and assume the model is right; since the same model made the cracks, recovering η shows consistency, not truth.
- **The rating:** every measured case is rated 5 "Excellent", both the 0.16 mm crack (the model says 64 % of the prestress is lost) and the 0.40 mm crack (79 % lost).
- **Serviceability** at η = 0.36 (case2 at GSD 0.02): truck-only midspan deflection 72.1 mm against 9.5 mm intact (Figure 3), band 69.1 – 72.1 – 74.4 mm. The limit L/800 is 30.6 mm: the intact girder uses 31 % of it, the damaged girder deflects 2.4 times the limit. The secant EI falls to 0.12 of the uncracked value at midspan, and 24 of 40 elements (x ≈ 4.9 – 19.6 m) are cracked by the passing truck (Figure 4).
- **Mesh:** the intact service deflection is 11.54 mm for N = 10 to 160 (spread 1e-9, round-off), as Hermite elements predict. With damage the relative difference to N = 160 is below 0.1 % for every N tested (0.0012 % from N = 80 to 160). N = 40 is used.
- **Branson cross-check** (same section and Mcr, only the interpolation changed): service deflection 73.3 against 81.5 mm with EC2 (Branson 10 % lower); truck-only 65.6 against 72.1 mm.
- **Monte Carlo stability:** the η percentiles for 1,000, 2,000 and 10,000 samples agree to about 0.001.
- **Sensitivity** (Figure 5): girder spacing 1.8 – 2.6 m moves η by 0.082, strand centroid 0.049, slab thickness 0.029, bond factor ξ 0.027, strand cover 0.014, girder fck 0.0009. The spacing alone moves η more than the measurement band at GSD 0.02 (0.046 wide): at a fine GSD the stand-in geometry, not the camera, limits the estimate.
- **Timing:** linear solve 1.5 ms, nonlinear analysis 34 ms, inverse 2.7 ms, Monte Carlo (2,000 samples) 5.2 s. A surrogate model is not needed.

![Figure 1. Predicted crack width at rest against the remaining prestress η (Stage-1 girder), with the Monte Carlo 5–95 % bands of the measured cases, the measured width ± u (grey) and the 3 px gate widths.](../results/figures/fig1_w_eta.png)

![Figure 2. All 12 synthetic cases: Monte Carlo band (blue), lower bound from a refusal at the 3 px gate (red), and refusals that give no bound at GSD 0.38 mm/px.](../results/figures/fig7_cases.png)

![Figure 3. Stage-1 girder: deflection under stage-B dead load + truck, intact (η = 1) and damaged (η = 0.36).](../results/figures/fig4_deflection.png)

![Figure 4. Stage-1 girder: secant EI / uncracked EI along the span under service load; the cracked zone emerges from each element's own moment.](../results/figures/fig3_ei_ratio.png)

![Figure 5. Sensitivity of the η estimate to the stand-in inputs (as-built 30-strand girder, case2 at GSD 0.02).](../results/figures/fig8_tornado.png)

## 6.3 Safety fixes in the pipeline link

- The parser of the pipeline's report now fails closed: a crack it cannot read stops the analysis instead of being skipped, and a crack the pipeline could not resolve blocks the governing width, because the unread crack may be the widest one (reporting a narrower crack would overestimate η, which is unconservative).
- A refusal is turned into a bound on η only when it is a "below resolution" refusal.
- The pre-cracked zone uses the moving-truck envelope, not one truck position.
- The figures were redrawn after review so that they show what was computed, not points that lie on the curve by construction.

**Numbers:** 256 tests passed and 3 skipped (the three that call the private pipeline); with the pipeline present they pass too. Stage 1 was merged on 12 October and tagged as the fallback deliverable.

# 7 Blocks 4–5: Stage 2, the grillage of the whole deck (15 October)

**Planned:** put the Stage-1 girders into a grillage of the S01 deck so that the FE model sits inside the secant loop and the inverse; invert η from a crack on the exterior girder G1 and on the interior girder G3; report the load shed to the neighbours and the deflection of every girder.

## 7.1 The model

- **Grid element** (Logan §5.4): 3 DOF per node (deflection w, torsional rotation θx, bending rotation θz); bending block from the Stage-1 beam element, torsion GJ/L, and the transformation k = Tᵀ·k̂·T for members in any direction in the plane.
- **Deck:** 5 girder lines at 2.20 m, 20 elements per girder; slab strips at every girder station with Hambly's properties (I = b·t³/12, J = b·t³/6); end and midspan diaphragms modelled as their webs (1.10 × 0.25 m); an overhang of 0.914 m.
- **Loads:** stage A on each bare girder from statics, as in Stage 1. Stage B on the grillage: dead load split to the girder lines by the lever rule; the HS20-44 × 1.30 truck as wheel lines, with one truck or two side by side in non-overlapping 12 ft design lanes (AASHTO Standard Specifications Art. 3.6 as quoted by the sources). Each girder's critical placement is held fixed for every η.
- **Nonlinear analysis:** every girder element follows the Stage-1 stiffness rule from its own moment; the cracked girder has η < 1, the others η = 1, and a neighbour that cracks under the shed load is handled by the same rule. Stage B is applied in 10 load steps, and new cracks are flagged one event at a time (largest M/Mcr first) from converged states only.
- **Inverse and Monte Carlo:** w_pred(η) is computed exactly on a 99-point η grid and joined linearly; the Monte Carlo runs on that table, and the exact inverse at the same widths confirms the percentiles to 1.4e-4 in η.

## 7.2 Verification

- The grid element reproduces Logan's Example 5.5 (d1y = −2.83 in, φ1x = 0.0295 rad, φ1z = −0.0169 rad and every printed member force) to Logan's three figures, and a symbolic derivation exactly.
- Identical girders under equal load each carry qL²/8 and deflect 5qL⁴/384EI; a single girder line reproduces the Stage-1 girder; a deck without torsion follows Courbon's rule with rigid cross members and the lever rule with flexible ones, and the real deck lies between them.
- An independent program, PyNite 2.0.2, solves the same linear intact grillage with a largest difference of 2.66e-11 mm in deflection (largest deflection 11.89 mm).

## 7.3 Results

- **Design:** strands G1–G5 = 30, 29, 27, 29, 30. The truck midspan moment is largest on G1 (1652.8 kN·m) with two trucks side by side.
- **Blind zone:** no crack at rest above η_dec = 0.536 (G1) and 0.615 (G3), against 0.557 for the Stage-1 girder.
- **η per measured case** (Monte Carlo 5 – 50 – 95 %). The synthetic cracks were made with the Stage-1 model, so the generating η is a 1-D value, not a truth for the grillage:

| Case | η used to generate (1-D) | Stage 1 | Grillage, G1 cracked | Grillage, G3 cracked |
|---|---|---|---|---|
| case2, GSD 0.02 | 0.36 | 0.338 – 0.361 – 0.385 | 0.263 – 0.287 – 0.314 | 0.330 – 0.353 – 0.378 |
| case2, GSD 0.05 | 0.36 | 0.312 – 0.367 – 0.433 | 0.236 – 0.294 – 0.374 | 0.304 – 0.359 – 0.434 |
| case3, GSD 0.02 | 0.21 | 0.196 – 0.215 – 0.234 | 0.120 – 0.138 – 0.157 | 0.192 – 0.209 – 0.227 |
| case3, GSD 0.05 | 0.21 | 0.169 – 0.214 – 0.262 | 0.093 – 0.138 – 0.185 | 0.166 – 0.209 – 0.255 |
| case3, GSD 0.10 | 0.21 | 0.124 – 0.214 – 0.312 | 0.050 – 0.137 – 0.236 | 0.123 – 0.208 – 0.304 |

- **The same crack means a lower η on G1.** As G1 softens it sheds part of its stage-B dead load to its neighbours: at η = 0.287 its stage-B dead-load midspan moment is 60.8 kN·m against 325.5 kN·m intact, so it needs less remaining prestress to open the same crack (Figure 6). At GSD 0.02 the choice of structural model moves η (1-D → grillage, G1: 0.074) more than the measurement does (band 0.051 wide). G3 stays close to Stage 1 (0.353 against 0.361), but it also has fewer strands (27 against 30), so this agreement is not a check of either model.
- **Load sharing** with G1 cracked at η = 0.287 (Figure 7): G1 loses 71.6 % of its stage-B service moment; G2 takes 59.6 % of it and G3 34.6 %. G2 cracks in 2 elements under the shed load, so its zero-tension design holds only for an intact deck. With G3 cracked at η = 0.353, G3 loses 88.8 % and G2 and G4 take 28.7 % each.
- **Deflection** (Figure 8): truck-only deflection of G1 9.52 → 21.9 mm at η = 0.287, inside L/800 = 30.6 mm. At the same η = 0.36 as Stage 1, G1 deflects 20.4 mm and G3 9.15 mm, against 72.1 mm for the Stage-1 girder. **The L/800 verdict reverses.** The Stage-1 girder cannot shed: its moment comes from statics with a fixed distribution factor (S/5.5), so it keeps its whole truck share however soft it becomes; in the grillage the softened girder passes the truck to its neighbours.
- **Model sensitivity** (G1, η = 0.287): the diaphragm stand-in moves the result more than the lane rule.

| Diaphragm, lane rule | G1 loses (%) | w_pred (mm) | G2 cracked elements | η median, case2 at GSD 0.02 |
|---|---|---|---|---|
| web, design lanes (as run) | 71.6 | 0.1604 | 2 | 0.287 |
| none, design lanes | 61.6 | 0.1948 | 1 | 0.313 |
| composite, design lanes | 76.0 | 0.1499 | 1 | 0.279 |
| web, clearance widths touching | 70.9 | 0.1604 | 2 | 0.287 |

  The verdict (inside L/800) holds for every variant: G1 truck deflection 20.5 – 23.9 mm.
- **The secant loop is now real:** every damaged analysis needs more than one pass (at least 5 per load step, up to 16 in the last service step of G1), with the moment change falling steadily (Figure 9). Proposal problem P5 is solved on the indeterminate model.
- **Torsion:** halving GJ of the cracked girder elements changes G1's w_pred by −0.07 % and its moment by +1.4 %.
- **Mesh** (Figure 10, relative to N = 80): w_pred at N = 20 within 1.5e-5 (G1) and 4.6e-5 (G3); G1's damaged service deflection within 1.2 % (not monotone). The inverse quantity is converged; G1's deflection carries about 1 % mesh uncertainty, which does not change the L/800 verdict.
- **Timing:** grillage build 4.5 s, linear solve 11.4 ms, nonlinear service analysis 2.37 s; the whole Stage-2 run takes 7.2 minutes.

![Figure 6. Predicted crack width at rest against η for the Stage-1 girder, the grillage with G1 cracked and with G3 cracked; dashed lines mark η_dec.](../results/figures/fig9_grillage_w_eta.png)

![Figure 7. G1 cracked at η = 0.287: stage-B service midspan moments of the five girders, intact and damaged.](../results/figures/fig11_moment_share_G1.png)

![Figure 8. G1 cracked at η = 0.287: truck-only midspan deflection of the five girders against L/800.](../results/figures/fig12_deflection_G1.png)

![Figure 9. Secant iteration on the grillage: relative moment change per stiffness update in the last load step.](../results/figures/fig13_secant_history.png)

![Figure 10. Grillage mesh study: relative change against N = 80.](../results/figures/fig14_grillage_mesh.png)

## 7.4 Changes forced by the grillage

- **Path-dependent cracking.** In the first version the set of cracked elements depended on the iteration path and on the load-step size. Stage B is now applied in increments and cracks are flagged one event at a time from converged states, so the answer no longer depends on how it was reached.
- **Lane rule.** The first truck placement let two trucks' clearance widths touch anywhere; it now follows the AASHTO design-lane rule, and the old rule is kept as a sensitivity.
- **Diaphragms.** A T-section diaphragm with the slab strip as flange made G1's result drift by 3 % per mesh doubling (the flange width is the element length); the web alone is mesh-independent and is bounded by the "none" and "composite" variants.
- **Labels and comparisons.** The η used to make the synthetic cracks is labelled as a 1-D value, never "η true", in every Stage-2 output, and the Stage-1/Stage-2 deflection comparison was redone at the same η, because the first version compared different η and different girders.

**Numbers:** 401 tests pass (398 and 3 skipped without the private pipeline). Stage 2 was merged on 15 October.

# 8 Block 6: Closing the proposal's outcomes (24 October)

**Planned:** the proposal compliance audit of 18 October (`docs/PROPOSAL_COMPLIANCE_AUDIT.md`) found three promised results that the report did not yet show: the response plotted against the measured width (Expected Outcome 2), a number for one cracked element against a damaged zone (Expected Outcome 3), and a hand calculation of one element inside the assembled system (Expected Outcome 1). It also asked for one table that maps every promise to its delivery.

**Changed and why:** a separate script, `scripts/run_proposal_checks.py`, computes the three results with the Stage-1 and Stage-2 models unchanged and writes `results/proposal_checks.json` and Figures 11 to 13; `stage1.json` and `stage2.json` are not touched. The Stage-1 girder got one new option, `damaged`, which lets only a chosen set of elements crack and keeps the others at their uncracked stiffness. Section 8.4 is the compliance table.

## 8.1 Response against the measured width (Expected Outcome 2)

The proposal asked for the midspan deflection and the stiffness ratio α = Ie/Ig plotted against the measured crack width, with the measurement band carried through. Each point of Figures 11 and 12 is one η: the model gives the crack width at rest w(η) and the response under the truck at the same η, so each curve is the response that a measured width implies under that model.

- **Deflection (Figure 11).** On the Stage-1 girder the truck-only midspan deflection is 9.5 mm intact, 38.0 mm at η_dec, 72.1 mm at w = 0.16 mm and 86.4 mm at w = 0.79 mm. The grillage girders stay far lower: G1 rises from 16.0 mm at its η_dec to 22.8 mm, G3 from 8.1 to 9.3 mm. Beyond w = 0.3 mm every curve rises by less than 7 %, so a wider crack adds little deflection.
- **The blind zone in deflection.** The dotted segment at w = 0 covers η = 1 to η_dec, where no crack is open at rest. On the Stage-1 girder the truck deflection already passes L/800 inside it (38.0 against 30.6 mm at η_dec): under this model a girder can fail the deflection limit before any crack is visible at rest. On the grillage it does not (G1 16.0 mm).
- **Stiffness ratio (Figure 12), two definitions.** (a) The midspan element ratio EI/EI_I under service load (stage-B dead load and truck), taken as the smaller of the two elements that meet at midspan: a local value, the counterpart of Branson's Ie/Ig at the cracked section. (b) The equivalent uniform ratio α = δ(η = 1)/δ(η): the factor by which a uniform EI would have to fall to give the same truck deflection, which is what α meant in the proposal's one-girder model. It is given for the single girder only, because on the grillage a deflection ratio also contains the load shed to the neighbours. At w = 0.16 mm the girder's midspan ratio is 0.122 and α = 0.132; at η_dec they are 0.208 and 0.250.
- **The proposal's limiting check.** A very wide crack should stop softening the girder. It does: above w = 0.4 mm every ratio stays between 0.100 and 0.113 (the girder's midspan ratio is 0.110 at 0.43 mm and 0.109 at 0.79 mm) and the deflection curves flatten. In the EC2 interpolation the cracked share ζ = 1 − β(Mcr/M)² approaches 1 as η and Mcr fall, so the response approaches that of the fully cracked section.
- **The measurement band.** The two Stage-1 cases at GSD 0.02 mm/px and the grillage reference case are drawn at w ± u with the Monte Carlo 5–95 % band of the response. At w = 0.16 ± 0.04 mm the girder deflects 69.1 – 72.1 – 74.4 mm, its α is 0.128 – 0.132 – 0.138 and its midspan ratio 0.119 – 0.122 – 0.126; G1 deflects 21.1 – 21.9 – 22.2 mm and G3 9.11 – 9.16 – 9.20 mm.

**The proposal's end points (added 1 November).** The proposal planned to solve again at the two ends of the measurement, w − u and w + u with u = 2·GSD. The table gives those two solves for every measured Stage-1 case beside the Monte Carlo band of Section 6.1. The end points always enclose the band and are about 20 % wider, because ±u is two standard deviations of the measurement (σ = u/2) while the 5–95 % band is ±1.645σ. The two methods agree; the band stays the main result because it also gives a median and carries through to the deflection and the stiffness ratio.

| Case | GSD (mm/px) | w ± u (mm) | η at w + u | η at w | η at w − u | Monte Carlo 5–95 % |
|---|---|---|---|---|---|---|
| case2 | 0.02 | 0.16 ± 0.04 | 0.334 | 0.361 | 0.390 | 0.338 – 0.385 |
| case2 | 0.05 | 0.15 ± 0.10 | 0.302 | 0.368 | 0.448 | 0.312 – 0.433 |
| case3 | 0.02 | 0.40 ± 0.04 | 0.193 | 0.215 | 0.237 | 0.196 – 0.234 |
| case3 | 0.05 | 0.40 ± 0.10 | 0.160 | 0.215 | 0.272 | 0.169 – 0.262 |
| case3 | 0.10 | 0.40 ± 0.20 | 0.111 | 0.215 | 0.334 | 0.124 – 0.312 |

![Figure 11. Truck-only midspan deflection against the crack width measured at rest, Stage-1 girder and grillage girders G1 and G3. Dotted at w = 0: η from 1 down to η_dec, where no crack is open at rest. Squares: measured cases at w ± u with the Monte Carlo 5–95 % band of the deflection.](../results/figures/fig15_width_deflection.png)

![Figure 12. Stiffness ratio against the crack width at rest: midspan element EI/EI_I for the three models and the equivalent uniform ratio α = δ(η = 1)/δ(η) for the Stage-1 girder, with the measured cases at GSD 0.02 mm/px.](../results/figures/fig16_width_stiffness.png)

## 8.2 One cracked element against a damaged zone (Expected Outcome 3)

The proposal put all the damage in one element at midspan and asked how much that element moves the global deflection compared with a fully iterated solve or a longer damaged zone. The Stage-1 girder answers it at η = 0.36 (the reference case) with three damage zones. Inside each zone the elements take their EC2 secant stiffness from their own moment, exactly as in the main model; all other elements keep their uncracked stiffness.

- **Tied to the mesh:** only the two elements that meet at midspan crack, a length of 2L/N. Every mesh here has a node at midspan, so the proposal's single element of length L/N at N = 10, 20 and 40 (2.45, 1.225 and 0.6125 m) has the same length as this zone at N = 20, 40 and 80.
- **Fixed physical length:** the central ℓ = L/5 = 4.90 m, which every mesh resolves exactly.
- **Distributed (the model of this report):** every element takes its stiffness from its own moment, so the cracked zone is wherever the moment exceeds the cracking moment. At η = 0.36 under service load it is 14.7 m long (x = 4.9 to 19.6 m).

| Truck-only midspan deflection (mm) | N = 10 | N = 20 | N = 40 | N = 80 | N = 160 |
|---|---|---|---|---|---|
| Two elements at midspan (ℓ = 2L/N) | 41.89 (ℓ 4.90 m) | 27.14 (2.45 m) | 18.70 (1.23 m) | 14.19 (0.61 m) | 11.87 (0.31 m) |
| Fixed ℓ = 4.90 m | 41.89 | 41.89 | 41.88 | 41.88 | 41.88 |
| Distributed | 72.17 | 72.13 | 72.12 | 72.12 | 72.12 |

The intact girder deflects 9.51 mm and L/800 is 30.6 mm.

**The number the proposal asked for.** The proposal's one-element model at its own meshes N = 10, 20 and 40 gives 27.1, 18.7 and 14.2 mm, against 72.1 mm for the distributed model. It captures 28 %, 15 % and 7.5 % of the extra deflection the damage causes (72.1 − 9.5 mm), about half as much at each mesh refinement, so the answer is set by the mesh and not by the crack. At N = 10 it would have put the damaged girder inside L/800 (27.1 against 30.6 mm) where the distributed model puts it at 2.4 times the limit.

**A longer zone.** With a fixed physical length the result no longer depends on the mesh (41.89 to 41.88 mm from N = 10 to 160). On the finest mesh the deflection grows with the damaged length (Figure 13, right) and reaches the distributed result exactly once the zone covers the 14.7 m that the moments crack; a longer zone adds nothing, because the elements outside it are uncracked anyway. The damaged length is therefore an output of the moment distribution, which is why the design review replaced the one-element model (P2 in Section 3).

![Figure 13. Stage-1 girder at η = 0.36. Left: truck-only midspan deflection against the number of elements for damage in the two elements at midspan (ℓ = 2L/N), in a fixed 4.90 m zone, and distributed by the moments. Right: N = 160, damage confined to a length ℓ around midspan.](../results/figures/fig17_damaged_length.png)

## 8.3 Hand calculation (Expected Outcome 1)

Appendix B works one element of the Stage-1 girder by hand with N = 10: the 4 × 4 stiffness matrix in numbers, the consistent load vector, the element-to-global DOF map, the assembled block at one node, the partition and the equilibrium check, each compared with the code. Under the stage-B dead load the FE midspan deflection is 2.0338 mm, equal to 5qL⁴/384EI, and each reaction is 52.95 kN, half the total load of 105.89 kN.

## 8.4 Proposal compliance

| Proposal item | Delivered | Where | Status |
|---|---|---|---|
| EO1: a verified 1-D beam FE solver, checked against a closed form and a hand calculation | Hermite beam element, work-equivalent loads, direct assembly, partition solve; nodal deflections exact for N = 1 – 160; one element worked by hand through assembly and equilibrium | 5.2, 5.3, 10, App. B | Delivered |
| EO2: midspan deflection and α = Ie/Ig against the measured width, with the measurement band | Figures 11 and 12 for the girder and the grillage; α as a midspan element ratio and as an equivalent uniform ratio; measured cases at w ± u with the Monte Carlo 5–95 % band | 8.1 | Delivered |
| EO3: one cracked element against a fully iterated solve or a longer damaged zone | truck deflection for damage at midspan (ℓ = 2L/N), in a fixed 4.90 m zone and distributed, N = 10 – 160, and against the damaged length | 8.2 | Delivered; the one-element result depends on the mesh |
| EO4: a clear line between the indicator and the DOH rating | the rating is reported beside η and the model never changes it | 1.2, 6.2 | Delivered |
| EO5, left for later: Branson iteration, as-built section, moving load, ML surrogate | Branson kept as a cross-check (EC2 is the main method, P1); a moving HS20-44 × 1.30 truck; the section is still a stand-in (P3); surrogate dropped because one solve takes milliseconds | 3, 4, 6 | Moving load done; as-built open; surrogate dropped |
| Method: N = 20 elements; partition instead of a penalty | 40 elements (girder), 20 per girder line (grillage); free block K_ff·D_f = F_f | 5.2, 7.1 | Delivered |
| Method: the crack width sets Ie of one element through Branson | width → strand stress increase (EC2 §7.3.4 inverted) → η → EC2 §7.4.3 stiffness in every element | 3 (P1, P2) | Replaced; reasons in P1 and P2 |
| Method: iterate Ie and Ma to convergence | one pass on the determinate girder (proved by a test); a real secant iteration on the grillage | 5.3, 7 | Delivered where the structure makes it meaningful (P5) |
| Method: solve again at w ± 2·GSD | η solved again at w − u and w + u (u = 2·GSD) for every measured case, beside a 2,000-sample Monte Carlo with σ = u/2 that gives 5–95 % bands of η, deflection and stiffness | 6.1, 6.2, 8.1 | Delivered; the band is the main result |
| Check: mesh N = 10, 20, 40 against the closed form | undamaged nodal deflections exact for every N; damaged girder N = 10 – 160, grillage N = 10 – 80 | 5.3, 6.2, 7.3 | Delivered and extended |
| Check: w → 0 gives the uncracked solution; a very wide crack stops softening | η = 1 gives the linear uncracked solution (test); w = 0 holds for every η from η_dec to 1, so a zero width does not identify an intact girder (blind zone); the ratios flatten at 0.100 – 0.113 above w = 0.4 mm | 6.2, 8.1 | Delivered; the blind zone is a finding |
| Case study: CPK72 S01, 24.5 m, illustrative section | span and supports from the DOH drawing; stand-in AASHTO Type IV composite girder and deck, every value tagged | 4 | Delivered as a stand-in case |

**Numbers:** 425 tests pass with the private pipeline configured (3 of them need it and are skipped without it). `run_proposal_checks.py` runs in about 3 minutes.

# 9 Where FEM sits in the project

| Code | Course method | Reference |
|---|---|---|
| `elements.beam_k`, `hermite` | beam element stiffness from Hermite shape functions | Logan (5th ed.) §4.1, Eqs. 4.1.7 and 4.1.14 |
| `elements.udl_equivalent`, `point_equivalent` | work-equivalent nodal loads {f0}; {f} = [k]{d} − {f0} | Ch. 3 slides 3-8 to 3-10; Logan §4.4, Table D-1 |
| `solve.solve`, `gridsolve.solve` | assembly, partition K_ff·D_f = F_f, reactions, member end forces | Ch. 3 slides 3-6 to 3-10; Logan §4.2–4.3, §5.4 |
| `elements.grid_k_local`, `grid_T`, `grid_k_global` | grid element: beam bending + torsion bar, k = Tᵀ·k̂·T | Logan §5.4, Eqs. 5.4.16–5.4.19, Example 5.5 |
| `nonlinear.solve_secant`, `secant_iterate` | secant stiffness algorithm, D(i+1) = Ks(D(i))⁻¹·F | Ch. 11 slides, "Nonlinear analysis algorithms", 3) secant stiffness |
| `rootfind.illinois` | secant method kept inside a bracket (outer inverse loop) | Ch. 11 (secant idea) |

The concrete physics meets the FE model in one place: the stiffness function hands each element's secant EI, computed from its moment, to a generic secant solver that knows nothing about concrete. The element, assembly and solve modules take EI and GJ per element and return displacements and member forces.

# 10 Verification summary

| What | Checked against | Result |
|---|---|---|
| Beam element, uncracked girder | closed form | nodal deflections exact for N = 1 – 160 |
| Secant loop on a determinate girder | statics | exactly one pass (P5) |
| EC2 crack width | European Concrete Platform Example 7.3 | wk = 0.1846 against 0.184 mm printed |
| Cracked prestressed section | Mast (1998) Example 1 | c = 17.262 in, σp = 171.96 ksi (printed 17.26, 171.97) |
| Staged section | two-stage transformed-section hand calculation | within 1e-5 |
| Inverse | forward model | η recovered to 1e-5 |
| Grid element | Logan Example 5.5 and a symbolic derivation | Logan's three figures; exact |
| Linear grillage | PyNite 2.0.2 | 2.66e-11 mm |
| Grillage limits | Courbon's rule, lever rule, qL²/8, 5qL⁴/384EI | the deck lies between the limits; exact for the identical-girder cases |
| Monte Carlo | 1k / 2k / 10k samples; exact inverse | stable to 0.001; table vs exact 1.4e-4 |
| Mesh | N = 10 – 160 (girder), 10 – 80 (grillage) | below 0.1 % (girder); w_pred 1e-5, G1 deflection 1.2 % (grillage) |
| Hand calculation, one element (Appendix B) | 12EI/L³ … 4EI/L, qL/2 and qL²/12 by hand; 5qL⁴/384EI; qL/2 | equal to the code; deflection and reactions to 1e-9 |
| Damage option of the Stage-1 girder | damage everywhere = distributed model; no damage = linear uncracked solution | to 1e-12 and 1e-9 |
| Whole code | 425 automated tests | all pass |

# 11 Discussion

- **What the indicator adds to the rating.** The rating is 5 "Excellent" for every crack this girder can show, yet the model reads the same cracks as a loss of 64 % to 79 % of the prestress, a girder that deflects up to 2.4 times its limit when analysed alone, and neighbours that take over most of its load when analysed as a deck. A width-based rating treats a crack of this size as minor, which suits reinforced concrete; for a prestressed girder the same width needs a mechanical reading.
- **Where the uncertainty is.** At GSD 0.02 mm/px the measurement band on η is about 0.05 wide. The girder spacing alone moves η by 0.08, and the choice between a single girder and the grillage by 0.07. With field data, the as-built drawings would narrow the estimate more than a better camera would.
- **What the camera must resolve.** At the pipeline's field GSD (0.38 mm/px) nothing this girder can show is measurable, and a girder above η_dec shows nothing at all. The indicator needs close-range images (GSD around 0.02 to 0.05 mm/px) and can only see losses above about 44 %.

A single girder with a fixed distribution factor cannot shed load, so it overstates the deflection of a damaged girder by a factor of about 3.5 and understates the demand on its neighbours. Only the indeterminate model, solved with the secant loop, shows both. This is the part of the answer that needs FEM.

On the AI side, the pipeline provides the measurement and its uncertainty; the FE model provides the mechanics. A second AI, Claude, wrote the code, with the student's engineering decisions recorded in the decision journal and the final student review tracked separately. The independent checks in Section 10 support the numerical results; authorship alone does not.

For a bridge owner, the present output should be read as a screening indicator rather than a capacity assessment. Its most actionable finding is the measurement requirement: the field GSD of 0.38 mm/px cannot resolve the crack widths predicted by this stand-in girder, whereas close-range images near 0.02–0.05 mm/px can support the inverse. Before acting on the estimated η, the owner would need the as-built geometry, the history and cause of the crack, and an independent structural assessment.

# 12 Limitations

- Every structural value is a stand-in until the as-built drawings of CPK72 are available: girder type, strands, concrete strengths, deck layout, overhang and diaphragms.
- The crack widths are synthetic. The masks are straight, constant-width bars, so they test the pipeline link and its resolution gates, not its measurement accuracy. The same model generated and inverted them (Stage 1), and the grillage bands cannot be scored against a truth.
- Not modelled: differential shrinkage of the deck, the timing of the prestress losses (all placed before casting), creep beyond what η absorbs, shear cracking and dynamics.
- Stage A comes from statics on the bare girder, so prestress lost after the deck was cast is not restrained by the deck; modelling that would push the G1 estimate of η lower still.
- The diaphragms are modelled as their webs, rigidly connected: the three-case sensitivity shows the spread (G1 loses 61.6 – 76.0 %) but does not bound it.
- The truck is placed for maximum midspan moment, which underestimates the maximum midspan live-load deflection by 1.9 % (measured in review, DECISIONS 2026-10-12); L/800 is compared with the truck-only deflection.
- The pipeline's ±u is read as a 95 % Gaussian band; the η band covers measurement noise only.
- A crack seen at rest may be a residual opening from a past overload (debris preventing closure); η would then be underestimated.
- The bond factor ξ = 0.6 for pretensioned strand is unverified in open sources.

# 13 Conclusions and next steps

1. Under the adopted crack and structural models, a crack width measured at rest can be mapped to estimates of remaining prestress, stiffness, truck deflection and deck load sharing, with a band from the measurement. The rule-engine rating does not provide this mechanical interpretation.
2. All five problems found in the proposal were resolved. In particular, the FE model now sits inside a real nonlinear loop: the secant iteration needs at least five passes per load step on the grillage, against exactly one on the single girder.
3. The deck model changes the engineering verdict: the single-girder model says the damaged girder fails L/800 by a factor of 2.4; the grillage says it passes, because the load goes to the neighbours, which then need checking themselves.
4. Putting all the damage in one element, as the proposal planned, gives an answer set by the mesh (27.1 mm at N = 10, 14.2 mm at N = 40, against 72.1 mm with the cracked zone set by the moments): the damaged length has to come from the moments, not from the element size.
5. The method is limited by what the camera can see (nothing at the field GSD of 0.38 mm/px; a blind zone up to 44 % loss) and by the stand-in structure, which moves η more than the measurement does at fine GSD.

**Next steps:** run the November 2026 field images through the same chain; replace the stand-ins with the BRB-CP as-built values; use the held-out check planned in the design (the predicted cracked-zone length against the crack extent seen in the images), which only becomes meaningful with real images.

# 14 Student decisions and review before submission

The student made three material engineering decisions recorded during the project: staged construction was retained, the HS20-44 load was multiplied by the DOH factor 1.30, and the whole-deck grillage became the main result while the single girder remained the fallback and comparison model. These choices changed the decompression threshold, service response and interpretation of the same measured crack.

The implementation was produced with AI assistance, so submission requires a separate oral and code review by the student. The review record is `report/STUDENT_REVIEW.md`. It covers the beam and grid elements, assembly and partition, the secant loop, the inverse, the uncertainty statement, and the distinction between a synthetic consistency check and field validation. The student declared the review complete on 1 November 2026; the review record states which answers were written down.

# References

Full source table with clauses, pages and links: `docs/sources.md` (S01–S30). Main references:

1. Logan, D. L. *A First Course in the Finite Element Method*, 5th ed.: §4.1–4.4 (beam element, equivalent loads), §5.4 (grid element, Example 5.5).
2. Course slides 01203541, Ch. 3 (beam, equivalent nodal loads) and Ch. 11 (nonlinear analysis algorithms).
3. JRC, *Eurocode 2: Background & Applications — Worked Examples*, EUR 26566 EN (2014): EC2 §7.3.4 crack-width equations and coefficients.
4. European Concrete Platform, *EC2 Worked Examples* and *Eurocode 2 Commentary* (2017): Example 7.3; strand modulus.
5. Sétra, *Eurocode 2 — Application to Concrete Highway Bridges* (2007): EC2 Eq. 7.10 and Table 6.2.
6. Mancini, *EC2: Concrete — Case Study*, JRC Eurocodes campaign (2011): EC2 §7.4.3 stiffness interpolation and β.
7. Mast, R. F. "Analysis of Cracked Prestressed Concrete Sections: A Practical Approach", *PCI Journal*, July–Aug 1998.
8. PCI, *Bridge Design Manual*, Appendix B (2011): AASHTO Type IV girder.
9. Hueste et al., TxDOT/TTI Report 0-4751-1 Vol. 2 (2006): AASHTO Standard Specifications design example (loads, impact, distribution, barrier, overhang).
10. Yang, Miao and Ghosn, NCHRP Report 776, Appendix D.1: grillage torsion constants after Hambly.
11. Roeder, Barth and Bergman, NCHRP Web Document 46 (2002): L/800 live-load deflection limit.
12. JICA, *The Second Mekong International Bridge Construction Project, Final Report* (2000): DOH standard HS20-44 × 1.30.
13. Chaini and Suparp, *Kasem Bundit Engineering Journal* 6(1), 2016: Thai I-girder bridges, AASHTO practice.
14. Bischoff, "A plea for unified deflection calculation of reinforced concrete flexural members", CSCE 2018: Branson's equation.
15. PyNiteFEA 2.0.2 (Python finite element library): independent check of the linear grillage.

# Appendix A: Reproducing the results

Python 3.11 or newer.

- `python -m pip install -r requirements.txt`
- `python -m pytest -q`: 425 passed with the configured private pipeline (captured 1 November 2026); without that configuration the 3 live-pipeline tests are skipped
- `python scripts/run_all.py`: Stage 1: `results/stage1.json` and figures 1–8 (about 75 s)
- `python scripts/run_stage2.py`: Stage 2: `results/stage2.json` and figures 9–14 (about 7 minutes)
- `python scripts/run_proposal_checks.py`: Block 6: `results/proposal_checks.json` and figures 15–17 (about 3 minutes; it reads stage1.json and stage2.json)
- `python scripts/build_demo_data.py`: the data of the browser demo, `docs/demo_data.js` (about 5 minutes); the demo itself is at https://palarch.github.io/crackfem/ and runs no FEM

The first two scripts start from the committed pipeline outputs in `data/claimbridge/`, so neither needs the private pipeline. The design document (`docs/superpowers/specs/`), the two implementation plans (`docs/superpowers/plans/`), the progress log, the source table and the decision journal are in `docs/`.

# Appendix B: Hand calculation of one beam element

Stage-1 girder, N = 10 elements, element length Le = 24.5 / 10 = 2.45 m. The composite girder is intact, EI = 9,969,881 kN·m², and the load is the stage-B dead load alone (wearing surface and barriers), q = 4.3222 kN/m acting downward. Units are kN and m; v is up positive and θ counter-clockwise positive, as in the course notes. The numbers are reproduced by `hand_calc()` in `scripts/run_proposal_checks.py`, which also asserts the checks of steps 4 and 5.

**Step 1: element stiffness matrix.** [k] = (EI/Le³)·[12, 6Le, −12, 6Le; 6Le, 4Le², −6Le, 2Le²; −12, −6Le, 12, −6Le; 6Le, 2Le², −6Le, 4Le²], with EI/Le³ = 9,969,881 / 14.7061 = 677,940.7 kN/m:

| [k] | v1 | θ1 | v2 | θ2 |
|---|---|---|---|---|
| v1 | 8,135,288 | 9,965,728 | −8,135,288 | 9,965,728 |
| θ1 | 9,965,728 | 16,277,356 | −9,965,728 | 8,138,678 |
| v2 | −8,135,288 | −9,965,728 | 8,135,288 | −9,965,728 |
| θ2 | 9,965,728 | 8,138,678 | −9,965,728 | 16,277,356 |

The entries are 12EI/Le³ (kN/m), 6EI/Le² (kN), 4EI/Le and 2EI/Le (kN·m). The matrix is symmetric, and a rigid-body translation (v1 = v2, θ1 = θ2 = 0) gives no force: 12 − 12 = 0 in the v rows and 6Le − 6Le = 0 in the θ rows.

**Step 2: consistent load vector.** For the downward load, {f} = [−qLe/2, −qLe²/12, −qLe/2, +qLe²/12] = [−5.2947 kN, −2.1620 kN·m, −5.2947 kN, +2.1620 kN·m]. These are the work-equivalent loads of the course notes, f = +{f0}.

**Step 3: DOF map and assembly.** Element e joins nodes e and e + 1, so its DOFs (v1, θ1, v2, θ2) go to global DOFs 2e − 1, 2e, 2e + 1 and 2e + 2, counting from 1: element 1 to DOFs 1–4, element 5 to DOFs 9–12 and element 6 to DOFs 11–14. The midspan node 6 owns DOFs 11 and 12. At node 2 the end of element 1 and the start of element 2 add:

| Node 2 block | v2 (DOF 3) | θ2 (DOF 4) |
|---|---|---|
| v2 (DOF 3) | 12EI/Le³ + 12EI/Le³ = 16,270,577 | −6EI/Le² + 6EI/Le² = 0 |
| θ2 (DOF 4) | 0 | 4EI/Le + 4EI/Le = 32,554,712 |

The global load at an interior node is −qLe = −10.589 kN with zero moment, because the two element end moments cancel. At the supports the load is −5.2947 kN, with −2.1620 kN·m at node 1 and +2.1620 kN·m at node 11.

**Step 4: partition and solve.** The supports restrain v at node 1 (DOF 1) and node 11 (DOF 21); the other 20 DOFs are free. The code solves K_ff·D_f = F_f on that 20 × 20 block and gets:

- midspan deflection v6 = −2.0338 mm, against the closed form 5qL⁴/384EI = 5 · 4.3222 · 24.5⁴ / (384 · 9,969,881) = 2.0338 mm (equal to 1e-9);
- end rotation θ1 = −2.6564e-4 rad, against −qL³/24EI = −2.6564e-4 rad.

Hermite elements with work-equivalent loads reproduce the exact nodal values for any N, which is why the proposal's mesh check on the uncracked beam could not test convergence (P2).

**Step 5: reactions and equilibrium.** R = K·D − F at the restrained DOFs gives R1 = R21 = 52.947 kN. Their sum, 105.894 kN, equals the total load qL = 4.3222 · 24.5 = 105.894 kN, and by symmetry each is qL/2.
