# CrackFEM submission

ClaimBridge is the background thesis work and the upstream measurement source: it supplies a crack width and a resolution band from an image-derived mask and surface scale. The work delivered to the FEM course begins at that hand-over. It contains the verified beam elements and assembly, the inverse estimate of remaining prestress, the nonlinear five-girder grillage, load redistribution and deflection results, uncertainty and sensitivity studies, tests, source records and the final report package. The present evidence uses synthetic masks and a stand-in bridge; it demonstrates the complete computational chain and its resolution limits, not field accuracy or structural capacity.

## Read these first

Quickest look, no install: <https://palarch.github.io/crackfem/> (browser demo; set a crack width and see η, stiffness and truck deflection; values interpolated between FEM solutions saved by `scripts/build_demo_data.py`).

1. `report/crackfem_report_FINAL.pdf` — main report.
2. `results/key_results.csv` — compact comparison of measured cases in the single-girder, exterior-girder and interior-girder models.
   `results/proposal_checks.json` with figures 15–17 — the proposal's Expected Outcomes 2 and 3 and the hand calculation (report Section 8, Appendix B); `docs/PROPOSAL_COMPLIANCE_AUDIT.md` ends with their status.
3. `results/figures/fig0_workflow.png` — boundary and flow from ClaimBridge to the FEM outputs.
4. `report/ClaimBridge_Background_and_FEM_Interface.pdf` — standalone background note for a supervisor who did not work on ClaimBridge; the matching `.docx` is editable.
5. `results/run_record.json` — result provenance and verification snapshot.
6. `report/STUDENT_REVIEW.md` — the student's review record and sign-off.

## Reproduce the numerical results

```text
python -m pip install -r requirements.txt
python -m pytest -q
python scripts/run_all.py
python scripts/run_stage2.py
python scripts/run_proposal_checks.py
python scripts/build_submission_assets.py
```

`run_all.py` regenerates Stage 1 and figures 1–8. `run_stage2.py` regenerates Stage 2 and figures 9–14. `run_proposal_checks.py` reads both result files and writes `results/proposal_checks.json` and figures 15–17. The first two use the committed parsed pipeline outputs and do not require the private ClaimBridge repository. Live ClaimBridge tests require the ignored per-machine pipeline configuration.

## Submission status

Final submission, 1 November 2026: 425 tests pass, report 26 pages with the headline figure on the cover, browser demo at https://palarch.github.io/crackfem/. The student declared the review complete on 1 November 2026 (`report/STUDENT_REVIEW.md`, and the `Reviewed by student` lines in `docs/PROGRESS_LOG.md`). No field-data or capacity claim should be added without independent evidence.
