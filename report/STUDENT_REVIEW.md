# Student review record

Complete this record in your own words before submitting the report. Short answers are sufficient when they show the mechanism clearly. Do not sign a block that you have not read.

## Identity and scope

- Student: Palid Archavachalee (6914500959)
- Project: Crack-width-informed prestress and stiffness indicator
- Review date: 1 November 2026
- Code/result commit reviewed: fb10d9b
- Status: The student declared the review complete on 1 November 2026. An oral review in chat on the same day covered Block 1 and the first Block 2 question; the answers to the questions below were not written down.

## Block 1 — sources and assumptions

- Which inputs are from the real CPK72 bridge, and which are stand-ins?
- Which stand-in has the largest demonstrated effect on η?
- One source or assumption that you would replace first with the as-built drawings:
- I reviewed `docs/sources.md` and `src/crackfem/params.py`: Yes / No

## Block 2 — Stage 1 FEM

- Why does the secant loop require one pass on the simply supported single girder?
- How are equivalent nodal loads and boundary conditions used in the partitioned solve?
- Why can a visible crack at rest indicate prestress loss under this model, and what alternative cause could make that inference wrong?
- I reviewed the Stage-1 reading path in `docs/PROGRESS_LOG.md`: Yes / No

## Block 3 — ClaimBridge and uncertainty

- Where does the course repository begin relative to ClaimBridge's image AI?
- Why can GSD 0.38 mm/px not support this indicator for the present stand-in girder?
- What uncertainty is included in the 5–95% η interval, and what uncertainty is omitted?
- I reviewed `claimbridge.py`, `uq.py` and `scripts/run_all.py`: Yes / No

## Block 4 — grillage

- Why does the grillage need several secant passes when the single girder needs one?
- Why does the same measured crack give different η estimates on G1, G3 and the Stage-1 girder?
- When G1 sheds moment to G2, what must be checked on the real bridge?
- I reviewed the Stage-2 reading path in `docs/PROGRESS_LOG.md`: Yes / No

## Block 6 — the proposal's outcomes (report Section 8, Appendix B)

- Redo Appendix B by hand: one element's [k] and {f}, then the node 2 block. Which number did you get differently at first, and why?
- Why does damage in one element give 27.1 mm at N = 10 but 14.2 mm at N = 40, while a fixed 4.90 m zone gives 41.9 mm at every N?
- The Stage-1 girder deflects 38.0 mm under the truck at η_dec, with no crack open at rest. What does that mean for an inspector who sees no crack?
- Section 11 (Discussion) ends with a bridge-owner paragraph that was written by the AI on 18 October, not by you. Keep it, rewrite it in your own words, or delete it: Keep (student, 1 November 2026)
- I reviewed `girder1d.py` (`damaged`), `scripts/run_proposal_checks.py` and report Section 8: Yes / No

## Engineering judgement

- The result that most changed my understanding of the problem:
- One modelling choice I disagree with or would change:
- What I would tell the bridge owner in two sentences:

## Sign-off

I understand that the current masks are synthetic, the structural geometry is mostly stand-in, and the output is a screening indicator rather than a capacity assessment.

- Student confirmation: Yes (Palid Archavachalee, 1 November 2026)
- Name and date: Palid Archavachalee, 1 November 2026
