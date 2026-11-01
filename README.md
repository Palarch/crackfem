# crackfem — crack width → remaining prestress → stiffness (FEM course project)

01203541 Finite Element Analysis of Structures, Kasetsart University, 2026 · Palid Archavachalee (6914500959)

A measured crack width on a prestressed bridge girder is turned into a mechanical statement:
how much prestress is left (η), how stiff the girder still is along the span, and how far it
deflects under the design truck — with an uncertainty band from the measurement. Crack widths
come from the author's thesis inspection pipeline (private) run on synthetic crack images.

Stage 1 (complete) models one interior girder of a simply supported span, built in two
construction stages: the bare girder carries the prestress and the weight of the girder and the
wet slab; the composite girder carries the wearing surface, the barriers and the design truck
through the beam finite-element model. Stage 2 (complete) puts the same girders into a grillage
of the whole five-girder deck (Logan's grid element, slab strips, diaphragms): the stiffness of
every girder element follows its own moment, so a cracked girder sheds load to its neighbours,
and the inverse for η runs through this indeterminate model for a crack on the exterior girder
G1 and on the interior girder G3.

The progress report is `report/crackfem_report_SUBMISSION_CANDIDATE.pdf` (with the matching
`.docx`); `SUBMISSION_README.md` lists what to read first.

## Quickstart

Python 3.11 or newer (developed and tested on 3.13).

    python -m pip install -r requirements.txt
    python -m pytest -q
    python scripts/run_all.py
    python scripts/run_stage2.py

`run_all.py` regenerates `results/stage1.json` and fig1–fig8 in `results/figures/`, and
`run_stage2.py` (about 7 minutes) regenerates `results/stage2.json` and fig9–fig14, both from the
committed pipeline outputs in `data/claimbridge/`; neither needs the private pipeline. The three
tests that call the pipeline are skipped unless a pipeline configuration is available (next
section), and `scripts/run_claimbridge.py` (re-measuring the synthetic masks) needs it too.

## Pipeline configuration (optional)

The thesis pipeline is private and this repository contains no detail of its code base: only
files cross the boundary (a PNG mask and a surface JSON go in, a text report comes out), and
the command that runs the pipeline is configured per machine. To run the live tests or
`scripts/run_claimbridge.py`:

1. Copy `pipeline.example.json` to `pipeline.local.json` (git ignores that file).
2. Replace the placeholders: `"cwd"` is the directory from which the pipeline's crack-width
   command is run, and `"argv"` is that command as a list of strings, up to but not including
   the `--mask` and `--surface` options (the bridge appends them).
3. Or keep the file elsewhere and point the environment variable `CRACK_PIPELINE_CONFIG` to it
   (it takes precedence over `pipeline.local.json`).

## What is where

- `src/crackfem/` — the package: section, crack width, stiffness, beam element, assembly and
  solve, secant iteration, inverse for η, Monte Carlo, the file bridge to the pipeline, figures;
  Stage 2: grid element and solve, deck layout and loads, the nonlinear grillage, the PyNite
  cross-check, the tabulated w(η)
- `tests/` — one test file per module, checked against independent evidence
- `scripts/` — synthetic cases, pipeline run, full results
- `pipeline.example.json` — placeholder for the per-machine pipeline configuration
- `data/claimbridge/` — synthetic crack masks and the parsed pipeline reports
- `results/` — generated numbers and figures

## Where to read

- Design and the five changes from the proposal: `docs/superpowers/specs/`
- Implementation plan: `docs/superpowers/plans/`
- What happened, block by block, with a reading guide and the results: `docs/PROGRESS_LOG.md`
- Every input value and its source: `docs/sources.md`
- Decisions and rejected alternatives: `docs/DECISIONS.md`

## Authorship

The code was written by Claude (Anthropic's AI) and is reviewed by the student; each entry in
`docs/PROGRESS_LOG.md` records what the student has reviewed and when.

## History

The project was developed in a local repository from 6 October to 1 November 2026. For
publication that history was squashed into one commit, because early commits held machine paths
and details of the private pipeline. Commit hashes quoted in `docs/` refer to the local history;
the dated entries in `docs/PROGRESS_LOG.md` and `docs/DECISIONS.md` record what each one did.

## Limits

All structural values are stand-ins until the as-built drawings are available (a standard
AASHTO Type IV girder with a composite slab, interior-girder load rules), and crack widths are
synthetic until the November 2026 field survey (see the spec §5.6). The synthetic masks test
the pipeline link and its resolution gates, not its measurement accuracy, and the same model
generates and inverts them. The synthetic widths were generated with the Stage-1 single-girder
model, so in Stage 2 the η they were made with is a 1-D model value, not a truth for the
grillage. The full list of limitations is in `docs/PROGRESS_LOG.md`.
