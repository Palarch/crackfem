# xray-ai pass on the report draft (2026-10-15)

Data record of one xray-ai pass: what the script found, what was changed, and what it found afterwards. Generated from `changes.json`, `xray_before.json` and `xray_after.json` in this folder.

| | Before | After |
|---|---|---|
| Script verdict | BLOCKED | BLOCKED |
| BLOCKER (placeholder-var) | 8 | 8 |
| TELL (em dash) | 39 | 3 |
| Pages | 17 | 17 |

Not changed on purpose:

- The 8 BLOCKERs are `{f0}`, `{f}`, `{d}`: the course's vector notation (Ch. 3 slides), not template variables. Kept.
- The 3 em dashes left are inside the titles of cited sources (references 3, 5, 6). Kept.
- The four highlighted STUDENT placeholders and the "Draft" label are for the student and were not touched.

Layout (not in the script): every figure paragraph is now kept with its caption (Word "Keep with next"); before, the Figure 3 caption fell onto the next page. Page renders: `before/pages_*.png`, `after_pages_*.png`.

## Changes

| # | Rule | Location | Before | After |
|---|---|---|---|---|
| 1 | em-dash (TELL) | heading 3 | # 3 Block 0 — Design review (6–7 October) | # 3 Block 0: Design review (6–7 October) |
| 2 | em-dash (TELL) | heading 4 | # 4 Block 1 — Source pack (7 October) | # 4 Block 1: Source pack (7 October) |
| 3 | em-dash (TELL) | heading 5 | # 5 Block 2 — Stage 1 core: one girder (7–12 October) | # 5 Block 2: Stage 1 core, one girder (7–12 October) |
| 4 | em-dash (TELL) | heading 6 | # 6 Block 3 — Link to the thesis pipeline | # 6 Block 3: Link to the thesis pipeline |
| 5 | em-dash (TELL) | heading 7 | # 7 Blocks 4–5 — Stage 2: the grillage of the whole deck | # 7 Blocks 4–5: Stage 2, the grillage of the whole deck |
| 6 | em-dash (TELL) | heading Appendix A | # Appendix A — Reproducing the results | # Appendix A: Reproducing the results |
| 7 | em-dash (TELL) | 5.1 item 2 | Stage A — prestress, girder and wet slab — acts on the bare girder; stage B — wearing surface, barriers and truck — acts on the composite girder. | Stage A (prestress, girder and wet slab) acts on the bare girder; stage B (wearing surface, barriers and truck) acts on the composite girder. |
| 8 | em-dash (TELL) | caption Figure 1 | ![Figure 1 —  | ![Figure 1.  |
| 9 | em-dash (TELL) | caption Figure 2 | ![Figure 2 —  | ![Figure 2.  |
| 10 | em-dash (TELL) | caption Figure 3 | ![Figure 3 —  | ![Figure 3.  |
| 11 | em-dash (TELL) | caption Figure 4 | ![Figure 4 —  | ![Figure 4.  |
| 12 | em-dash (TELL) | caption Figure 5 | ![Figure 5 —  | ![Figure 5.  |
| 13 | em-dash (TELL) | caption Figure 6 | ![Figure 6 —  | ![Figure 6.  |
| 14 | em-dash (TELL) | caption Figure 7 | ![Figure 7 —  | ![Figure 7.  |
| 15 | em-dash (TELL) | caption Figure 8 | ![Figure 8 —  | ![Figure 8.  |
| 16 | em-dash (TELL) | caption Figure 9 | ![Figure 9 —  | ![Figure 9.  |
| 17 | em-dash (TELL) | caption Figure 10 | ![Figure 10 —  | ![Figure 10.  |
| 18 | em-dash (TELL) | ref 1 | 5th ed. — §4.1 | 5th ed.: §4.1 |
| 19 | em-dash (TELL) | ref 3 | EUR 26566 EN (2014) — EC2 | EUR 26566 EN (2014): EC2 |
| 20 | em-dash (TELL) | ref 4 | Commentary* (2017) — Example | Commentary* (2017): Example |
| 21 | em-dash (TELL) | ref 5 | Highway Bridges* (2007) — EC2 | Highway Bridges* (2007): EC2 |
| 22 | em-dash (TELL) | ref 6 | campaign (2011) — EC2 | campaign (2011): EC2 |
| 23 | em-dash (TELL) | ref 8 | Appendix B (2011) — AASHTO | Appendix B (2011): AASHTO |
| 24 | em-dash (TELL) | ref 9 | Vol. 2 (2006) — AASHTO | Vol. 2 (2006): AASHTO |
| 25 | em-dash (TELL) | ref 10 | Appendix D.1 — grillage | Appendix D.1: grillage |
| 26 | em-dash (TELL) | ref 11 | Document 46 (2002) — L/800 | Document 46 (2002): L/800 |
| 27 | em-dash (TELL) | ref 12 | Final Report* (2000) — DOH | Final Report* (2000): DOH |
| 28 | em-dash (TELL) | ref 13 | 6(1), 2016 — Thai | 6(1), 2016: Thai |
| 29 | em-dash (TELL) | ref 14 | CSCE 2018 — Branson | CSCE 2018: Branson |
| 30 | em-dash (TELL) | ref 15 | library) — independent | library): independent |
| 31 | em-dash (TELL) | Appendix A, pytest | `python -m pytest -q` — 401 | `python -m pytest -q`: 401 |
| 32 | em-dash (TELL) | Appendix A, run_all | `python scripts/run_all.py` — Stage 1 | `python scripts/run_all.py`: Stage 1 |
| 33 | em-dash (TELL) | Appendix A, run_stage2 | `python scripts/run_stage2.py` — Stage 2 | `python scripts/run_stage2.py`: Stage 2 |
| 34 | judge: repeated template heading | 5.4 | ## 5.4 Problems met and fixed | ## 5.4 What the reviews caught |
| 35 | judge: repeated template heading | 6.3 | ## 6.3 Problems met and fixed | ## 6.3 Safety fixes in the pipeline link |
| 36 | judge: repeated template heading | 7.4 | ## 7.4 Problems met and fixed | ## 7.4 Changes forced by the grillage |
| 37 | judge: machine-even bullets | 10, 4th bullet | - **Why the deck model matters.** A single girder | A single girder |
| 38 | judge: machine-even bullets | 10, 5th bullet | - **AI in the project.** The AI pipeline | On the AI side, the pipeline |
| 39 | judge: adjective instead of number | 11, 6th bullet | which slightly underestimates the maximum deflection; | which underestimates the maximum midspan live-load deflection by 1.9 % (measured in review, DECISIONS 2026-10-12); |

Files: `before/` holds the draft as it was (md, docx, PDF, page renders); the current draft is `report/report_draft.md` and `report/crackfem_report_DRAFT.docx`.
