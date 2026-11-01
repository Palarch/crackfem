# Proposal compliance audit

**Proposal reviewed:** *FEM–AI Structural Indicator for Crack-Width-Informed Deflection: A Case Study on Bridge CPK72 (Route 304)*, dated 15 September 2026  
**Project reviewed:** current workspace and submission candidate, 18 October 2026  
**Verification run:** `python -m pytest -q --basetemp=tmp/pytest` — **401 passed**

## Overall finding

The project meets the proposal's central purpose and contains substantially more FEM work than the minimum promised. It provides a verified one-dimensional Euler–Bernoulli girder, propagates a ClaimBridge crack-width measurement into structural quantities, carries measurement uncertainty, keeps the FEM indicator separate from the DOH rating, and adds a nonlinear five-girder grillage that demonstrates load redistribution.

It is not a literal implementation of every proposed assumption. The final model changed several methods after the design review found that the original width-to-Branson link and one-element damage model were not mechanically defensible. These changes are documented and improve the work, but the final submission should make the replacement of each promised item easy for the lecturer to audit.

## Compliance matrix

| Proposal commitment | Evidence in the current project | Status |
|---|---|---|
| ClaimBridge supplies crack width and its measurement band | Parsed pipeline outputs for 12 synthetic mask/GSD cases; file-only adapter; `claimbridge.py`; Sections 1.2 and 6 | Met for the software interface; field measurement accuracy is not validated |
| Verified 1-D Euler–Bernoulli FE solver | Hermite beam element, consistent nodal loads, direct assembly, partitioned solve, reactions and member forces; Sections 5, 8 and 9 | Met |
| Switch boundary conditions from a penalty method to partitioning | `solve.solve` and `gridsolve.solve` solve the free block `K_ff D_f = F_f` | Met |
| Check against closed-form beam response | Exact nodal deflections for the simply supported uniform-load case from N = 1 to 160; other statics and virtual-work checks | Met |
| Hand calculation of one element stiffness/load vector | Element stiffness and load vectors are checked symbolically in automated tests; a transformed-section hand calculation and Logan grid example are reported | Partly met: the promised numerical beam-element hand calculation is not shown in the report |
| Mesh study at N = 10, 20 and 40 | Stage 1 uses N = 10–160 and Stage 2 uses N = 10–80, including the damaged response | Met and extended |
| Crack width controls stiffness through Branson effective inertia | The design review shows that Branson has no crack-width input. The main chain now uses EC2 crack width and EC2 tension-stiffening/secant stiffness; Branson remains a quantitative cross-check | Replaced with a justified method |
| Iterate stiffness and response to convergence | A simply supported girder correctly requires one pass because its moments are statically determinate. The indeterminate grillage requires multiple secant passes and redistributes load | Met in Stage 2; the proposal's intended nonlinear mechanism is demonstrated there |
| One damaged element at midspan | Stiffness is assigned to every element from its own moment, producing a mesh-independent cracked zone | Replaced with a justified distributed-zone model |
| Quantify one-element damage versus a fully iterated or longer damaged zone | The report compares the single-girder and grillage responses and reports the emerged cracked-zone length, but it does not give the specific one-element-versus-physical-damaged-length comparison promised in Expected Outcome 3 | Outstanding or must be explicitly marked as superseded |
| Plot midspan deflection and stiffness ratio `alpha = Ie/Ig` against measured crack width | The report contains width-versus-remaining-prestress, spatial EI ratio, and deflected-shape figures. It does not contain direct deflection-versus-width and EI-ratio-versus-width curves | Outstanding |
| Propagate `w ± 2 GSD` into a bounded structural result | Measurement uncertainty is propagated by 2,000-sample truncated-normal Monte Carlo with `sigma = u/2`, and 5–95% intervals are reported | Met in purpose, changed statistically; deterministic endpoint results at `w-u` and `w+u` are not tabulated |
| Keep the structural indicator separate from the DOH rating | Report, code and ClaimBridge background note state that the rule-engine rating is retained and never changed by FEM | Met |
| Use CPK72 span S01, 24.5 m | Span and support condition are tied to the DOH drawing | Met |
| Use illustrative section/material/load values until as-built drawings are checked | The rectangle was replaced by a stand-in AASHTO Type IV composite girder and stand-in deck; all such inputs are labelled | Met as an explicitly provisional case study |
| Moving load and surrogate ML remain later extensions | A moving HS20-44 × 1.30 truck was implemented. The surrogate was dropped because measured solve times do not justify it | Moving load exceeded the proposal; surrogate omission is justified |
| Advisory output rather than a capacity or DOH-rating decision | The report repeatedly limits the result to a screening indicator based on synthetic masks and a stand-in structure | Met |

## Items to close before submission

### Priority 1 — required for a defensible hand-in

1. **Complete the student review record.** `report/STUDENT_REVIEW.md` is blank and seven entries in `docs/PROGRESS_LOG.md` remain marked pending. This is especially important because the report states that AI wrote the implementation.
2. **Add a one-page proposal compliance table to the report.** The design-review table explains why the method changed, but the lecturer should not have to infer whether every Expected Outcome was delivered.
3. **Add the two direct curves promised in Expected Outcome 2:** measured crack width versus truck midspan deflection, and measured crack width versus a clearly defined stiffness ratio. State whether the stiffness ratio is the midspan element value, the minimum value, or an equivalent global value.
4. **Resolve Expected Outcome 3 explicitly.** Either add a mesh-independent damaged-length sensitivity, such as fixed physical lengths around midspan, or state that the one-element comparison was superseded because its damaged length equals `L/N` and is therefore not a valid physical comparison. A numerical replacement should then be named, such as distributed-zone girder versus grillage load redistribution.
5. **Show the promised hand calculation.** Add one numerical beam-element stiffness matrix, one consistent load vector, the DOF mapping into the global system, and one reaction/equilibrium check. The existing symbolic tests support this appendix but do not replace a visible hand-worked example for assessment.

### Priority 2 — improves traceability

6. Add deterministic endpoint values from `w-u` and `w+u` beside the Monte Carlo interval, or explicitly record that the proposal's endpoint band was replaced by a 5–95% probabilistic interval and why.
7. Regenerate `results/run_record.json` after the final commit and final document build. Its current record says the working tree was dirty and predates the latest ClaimBridge background additions.
8. Keep the field-data and as-built limitations prominent. The current synthetic masks verify the software chain, not ClaimBridge measurement accuracy, and most bridge properties remain stand-ins.

## Submission judgement

The project is technically beyond the proposed minimum and the main FEM objective is achieved. It should be described as **proposal purpose achieved with documented engineering revisions**, rather than as an unchanged execution of the original method. The remaining weaknesses are mostly assessment evidence and presentation, with the direct width-response plots, the explicit treatment of Expected Outcome 3, the visible hand calculation, and the student's own review being the clearest gaps.

## Status update — 24 October 2026

The audit text above is kept as written on 18 October. Items closed on 24 October (report Section 8, Appendix B; DECISIONS 2026-10-24):

| Item | Status |
|---|---|
| 1. Student review record | **Open**: needs the student; `report/STUDENT_REVIEW.md` now also has a Block 6 section |
| 2. Compliance table in the report | Closed: report Section 8.4 |
| 3. Width → deflection and width → stiffness-ratio curves | Closed: Figures 11 and 12 (Section 8.1); the ratio is defined as the midspan element value and as an equivalent uniform value |
| 4. Expected Outcome 3 | Closed: Section 8.2, one element vs fixed length vs distributed, N = 10–160, and against the damaged length |
| 5. Visible hand calculation | Closed: Appendix B, asserted in `scripts/run_proposal_checks.py` |
| 6. Deterministic w − u / w + u endpoints | Open (Priority 2); the report states that the Monte Carlo band replaces the two end points |
| 7. Regenerate `run_record.json` after the final commit | Open (Priority 2); do it after the student review |
| 8. Keep limitations prominent | Unchanged |

## Status update — 1 November 2026

| Item | Status |
|---|---|
| 1. Student review record | **Open**: the student chose to finish the build first and do the review last |
| 6. Deterministic w − u / w + u endpoints | Closed: η solved again at w − u and w + u for all five measured Stage-1 cases (`end_point_rows` in `scripts/run_proposal_checks.py`, report Section 8.1 table); the end points enclose the Monte Carlo 5–95 % band and are about 20 % wider |
| 7. Regenerate `run_record.json` after the final commit | Open until the publication commit exists; then recorded in a commit of its own (a file cannot hold the hash of the commit that contains it) |

## Status update — 1 November 2026, final

| Item | Status |
|---|---|
| 1. Student review record | **Closed by the student's declaration**: the student declared the review complete on 1 November 2026 and asked for no further questions; `report/STUDENT_REVIEW.md` records the declaration, the sign-off and the Discussion decision (keep), and states that the written answers were not recorded. |
| 7. Regenerate `run_record.json` after the final commit | Closed: recorded for the final commit in a commit of its own. |
