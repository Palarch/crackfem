# Plan 1 — Stage 1: Crack-Width → Prestress → Stiffness Pipeline on One Girder

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and verify the complete Stage-1 chain on one simply supported prestressed girder (CPK72 S01 stand-in): synthetic crack mask → ClaimBridge crack-width command → inverse for remaining prestress η → secant-stiffness beam FE → deflection/stiffness indicator with a Monte Carlo band. This is the fallback deliverable that must exist by 2026-10-16.

**Architecture:** A small Python package `crackfem` (src layout) whose FE core (`elements`, `model`, `solve`, `nonlinear`) knows nothing about concrete; concrete physics lives in `section`, `crackwidth`, `stiffness`, `prestress`; `girder1d` wires them; `inverse`/`rootfind`/`uq` answer the question; `claimbridge` exchanges only files with the private thesis pipeline. Every numeric input lives in `params.py` with a source id that must exist in `docs/sources.md`.

**Tech Stack:** Python 3.13, numpy 2.3.5, scipy 1.17.1, sympy 1.14.0, pillow 12.0.0, pytest 8.3.5, matplotlib (installed in Task 0 and pinned). The ClaimBridge pipeline (local, read-only), run through the command configured per machine in `pipeline.local.json` (see `pipeline.example.json`).

**Spec:** `docs/superpowers/specs/2026-10-07-crack-informed-grillage-design.md`. **Plan 2** (B4 grillage, B5 grillage inverse/MC/tornado, B6 report) is written on 2026-10-19/20 from the B1 values and Stage-1 results; see the end of this file for its scope.

---

## As-built notes (added 2026-10-12, Plan 1 close-out)

Plan 1 was executed task by task with a two-stage review after each group. The code blocks below are the plan as written; the repository differs where reviews or student decisions changed it. Every change is a dated line in `docs/DECISIONS.md`; the main ones:

- **Student decisions (2026-10-12):** staged construction (Task 15b added; Tasks 9, 11, 12, 16, 20, 21 amended in this file) and HS20-44 × 1.30.
- **Plan test errors found by implementers and corrected:** Task 7 compared cracked/uncracked states outside the window where they can agree (now 300 kN·m); Task 13 used a round-off bound below the noise floor (now 1e-9) and summed a moment reaction into a force check (now vertical reactions only).
- **Hardening after reviews (code differs from the blocks below):** truck scan includes the exact kink placements; secant iteration also checks the stiffness residual; inverse uses bisection for the decompression edge, a residual check and NaN guards; Monte Carlo uses inverse-CDF sampling and reports dropped statuses; pre-cracked zone uses the moving-truck envelope; the ClaimBridge parser fails closed; a third scenario (η_dec − 0.35); refusal bounds (nominal and 95 %); as-built tornado (`Girder1D(n_strands=...)`); figures made honest (fig1, fig2, fig6, fig7, fig8); extra input validation in model/solve/crackwidth/claimbridge.
- **Convention correction:** `f_eq = +{f0}` (see Conventions).
- The README, `docs/PROGRESS_LOG.md` (Entries 0–3) and `results/stage1.json` describe the as-built state.


## Conventions (read once; every task relies on them)

- **Units:** SI base units everywhere inside the code — N, m, Pa, N·m. Convert only for display (`*_mm`, `*_kNm` keys).
- **Section:** `y` = height above the girder soffit (m). Tension positive. Curvature `kappa > 0` = sagging. Strain `eps(y) = eps0 − kappa·y` (so `eps0` is the soffit strain). Strand stress `sigma_p = Ep·(eps_p0 + eps(y_p))`, where `eps_p0` is the strand strain locked in when the surrounding concrete has zero strain. η scales `eps_p0` (η = 1 → as designed).
- **Beam FE:** DOFs per node `(v, θ)` → global indices `(2i, 2i+1)`. `v` up positive, `θ` counter-clockwise positive, loads up positive (gravity is negative). System `K D = F` with `F` = assembled consistent (equivalent) nodal loads. Element end forces `f_e = k_e d_e − f_eq,e` = `[V1, M1, V2, M2]` acting on the element. Internal moment (sagging +) at local `x`: `M(x) = −M1 + V1·x + q·x²/2 + Σ P·(x − a)` for point loads at `a < x`. The Ch3 slides' `{f0}` is the work-equivalent load itself, so `f_eq = +{f0}` and `{f} = [k]{d} − {f0}`; fixed-end reactions are `−{f0}` (corrected in Task 22 after checking the slides; the first draft of this plan said `f_eq = −F0`).
- **FE carries external loads only**, measured from a reference state: secant `EI = (M − M_ref) / (kappa(M) − kappa(M_ref))`. Camber is reported separately.
- **Staged construction (student decision 2026-10-12):** stage A (prestress + girder + wet slab) acts on the bare girder and its moments come from statics; stage B (wearing surface, barriers, truck) acts on the composite girder through the FE model, so `M_ref` = the stage-A moment. Each staged section carries the slab's locked-in strain `eps_free` (the stage-A strain plane at slab level). η scales today's strand strain; the slab strain locked in at construction (η = 1) does not change.
- **Commands** run from `fem_project/`. Tests: `python -m pytest ...`. Commits end with the trailer line `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Public-repo hygiene:** never commit ClaimBridge code, field data, absolute local paths, URLs of internal services, or keys. Only synthetic masks/surfaces and parsed pipeline output.

## File structure

| File | Responsibility |
|---|---|
| `pyproject.toml`, `requirements.txt`, `.gitignore`, `README.md` | packaging/test config, pinned deps, quickstart |
| `src/crackfem/params.py` | every numeric input (`Inputs`) + `SOURCES` map to `docs/sources.md` |
| `src/crackfem/materials.py` | EC2 fcm/fctm/Ecm; ACI Ec/fr (cross-check only) |
| `src/crackfem/geometry.py` | trapezoid outline, exact gross properties, fibre discretisation, area below a cut |
| `src/crackfem/section.py` | fibre section equilibrium (uncracked / cracked), Δσp, uncracked EI, Mcr |
| `src/crackfem/prestress.py` | zero-tension prestress sizing, locked-in strain |
| `src/crackfem/staging.py` | staged construction: bare-girder section, slab locked-in strain, staged prestress design |
| `src/crackfem/crackwidth.py` | EC2 §7.3.4 crack width and its closed-form inverse; hc,ef |
| `src/crackfem/stiffness.py` | EC2 §7.4.3 secant EI; ACI Branson secant EI (cross-check) |
| `src/crackfem/elements.py` | beam element k, Hermite shape functions, consistent loads |
| `src/crackfem/model.py` | beam model data (nodes, loads, supports) |
| `src/crackfem/solve.py` | assembly, partition solve, end forces, internal moments |
| `src/crackfem/loads.py` | dead load per girder, HS20-44 truck, impact, distribution fraction |
| `src/crackfem/nonlinear.py` | direct iteration with secant stiffness (Ch11) |
| `src/crackfem/rootfind.py` | Illinois bracketed secant |
| `src/crackfem/inverse.py` | η from crack width; decompression bound; monotonicity |
| `src/crackfem/girder1d.py` | Stage-1 wiring: w_pred(η), stiffness, indicator |
| `src/crackfem/uq.py` | Monte Carlo on w; tornado |
| `src/crackfem/claimbridge.py` | masks, surface JSON, run/parse the pipeline CLI |
| `src/crackfem/figures.py` | all figures from arrays |
| `scripts/make_synthetic_cases.py`, `scripts/run_claimbridge.py`, `scripts/run_all.py` | data generation, thesis-pipeline run, full results |
| `tests/test_*.py`, `tests/fixtures/*` | one test file per module; published values; real pipeline output |
| `docs/sources.md`, `docs/PROGRESS_LOG.md`, `docs/DECISIONS.md` | evidence trail and progress narrative |

---

## Task 0: Project scaffold, matplotlib, spec fix, progress log

**Files:**
- Create: `pyproject.toml`, `requirements.txt`, `.gitignore`, `src/crackfem/__init__.py`, `tests/__init__.py` (empty), `docs/PROGRESS_LOG.md`
- Modify: `docs/superpowers/specs/2026-10-07-crack-informed-grillage-design.md` (§5.1 concrete row), `docs/DECISIONS.md`

- [ ] **Step 1: Install matplotlib (approved dependency, 2026-10-07)**

Run: `python -m pip install "matplotlib>=3.10"` then `python -m pip show matplotlib`
Expected: `Name: matplotlib` and a `Version:` line. Use that exact version in Step 2.

- [ ] **Step 2: Write `requirements.txt`** (replace `X.Y.Z` with the version printed in Step 1)

```text
numpy==2.3.5
scipy==1.17.1
sympy==1.14.0
pillow==12.0.0
pytest==8.3.5
matplotlib==X.Y.Z
```

- [ ] **Step 3: Write `pyproject.toml`**

```toml
[project]
name = "crackfem"
version = "0.1.0"
description = "Crack-width-informed prestress and stiffness indicator (FEM course project)"
requires-python = ">=3.11"

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
```

- [ ] **Step 4: Write `.gitignore`**

```text
__pycache__/
*.pyc
.pytest_cache/
```

- [ ] **Step 5: Write `src/crackfem/__init__.py`**

```python
"""crackfem: crack-width-informed prestress and stiffness indicator (FEM course project)."""
```

and an empty `tests/__init__.py`.

- [ ] **Step 6: Fix spec §5.1 concrete row to the EC2 family**

In the spec table row `| Concrete | girder f'c 40 MPa, slab 28 MPa; E and fr from code formulas | ACI 318-19 Ch. 19 |` replace with:

```markdown
| Concrete | girder fck 40 MPa, slab 28 MPa; Ecm and fctm from EC2 Table 3.1 (main path); ACI Ec and fr only inside the Branson cross-check | EN 1992-1-1 Table 3.1; ACI 318-19 §19.2 |
```

Append to `docs/DECISIONS.md` under a new heading `## 2026-10-08 — Plan 1`:

```markdown
- **Concrete modulus and tensile strength from EC2 Table 3.1 (Ecm, fctm) on the main path**, because crack width and stiffness are EC2; ACI Ec and fr are used only in the Branson cross-check. Rejected: ACI values on the main path (mixes code families, spec §5.1 first draft).
- **Spec split into Plan 1 (B1–B3, Stage 1) and Plan 2 (B4–B6, grillage + report)**: Plan 2 is written after the B1 source pack and Stage-1 results so grillage numbers are not locked before they are verified.
```

- [ ] **Step 7: Write `docs/PROGRESS_LOG.md`**

```markdown
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

**Changed & why:** inverse for remaining prestress η from the crack width; EC2 crack width + stiffness; per-element stiffness from each element's own moment; staged 1-D girder → grillage so FE sits inside the nonlinear and inverse loops (spec §2–§5, DECISIONS 2026-10-07).

**Numbers:** none yet. A probe on 2026-10-07 confirmed the pipeline's crack-width command measures a synthetic 6 px crack at GSD 0.05 mm/px as 0.30 ± 0.10 mm and the rule engine rates it 5 (Excellent).

**Reading guide:** spec §1–§5.

**Reviewed by student:** spec approved in chat 2026-10-07.
```

- [ ] **Step 8: Verify pytest runs (no tests yet)**

Run: `python -m pytest`
Expected: `no tests ran` (exit code 5 is fine at this point).

- [ ] **Step 9: Commit**

```bash
git add pyproject.toml requirements.txt .gitignore src tests docs
git commit -m "chore: scaffold crackfem package, pin deps, progress log entry 0" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 1: `params.py` + `docs/sources.md` skeleton

**Files:**
- Create: `src/crackfem/params.py`, `docs/sources.md`
- Test: `tests/test_params.py`

- [ ] **Step 1: Write the failing test** — `tests/test_params.py`

```python
import re
from pathlib import Path

from crackfem.params import SOURCES, Inputs, field_names

SOURCES_MD = Path(__file__).resolve().parents[1] / "docs" / "sources.md"


def _flatten(value):
    if isinstance(value, tuple):
        for item in value:
            yield from _flatten(item)
    else:
        yield value


def test_every_input_has_a_source_id():
    assert field_names() == set(SOURCES)


def test_every_source_id_is_a_row_of_sources_md():
    text = SOURCES_MD.read_text(encoding="utf-8")
    rows = set(re.findall(r"^\|\s*(S\d\d)\s*\|", text, flags=re.MULTILINE))
    assert set(SOURCES.values()) <= rows


def test_all_numeric_inputs_are_positive():
    inp = Inputs()
    for name in field_names():
        assert all(v > 0 for v in _flatten(getattr(inp, name))), name
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_params.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'crackfem.params'`

- [ ] **Step 3: Write `src/crackfem/params.py`**

```python
"""Every numeric input of the model, in SI base units (N, m, Pa).

Each field of Inputs has a source id in SOURCES; every id is a row of
docs/sources.md. Rows tagged [stand-in] are not the real CPK72 girder.
Values are candidates until Task 2 (source pack) verifies them.
"""
from dataclasses import dataclass, fields

KIP = 4448.2216152605       # N per kip
FT = 0.3048                 # m per ft
INCH = 0.0254               # m per inch
KSI = 6.894757293168361e6   # Pa per ksi


@dataclass(frozen=True)
class Inputs:
    # span and deck
    span: float = 24.50
    girder_spacing: float = 2.20
    slab_thickness: float = 0.20
    wearing_thickness: float = 0.05
    n_girders: int = 5
    n_barriers: int = 2
    barrier_weight: float = 7.0e3
    # girder outline, bottom to top: (height_in, width_bottom_in, width_top_in)
    girder_outline_in: tuple = ((8.0, 26.0, 26.0), (9.0, 26.0, 8.0), (23.0, 8.0, 8.0),
                                (6.0, 8.0, 20.0), (8.0, 20.0, 20.0))
    # strands
    strand_diameter: float = 12.7e-3
    strand_area: float = 98.7e-6
    fpu: float = 1860e6
    Ep: float = 195e9
    fpj_ratio: float = 0.75
    lump_sum_loss: float = 45.0 * KSI
    strand_centroid: float = 0.10
    strand_cover: float = 0.045
    # concrete
    fck_girder: float = 40e6
    fck_slab: float = 28e6
    unit_weight_concrete: float = 24.0e3
    unit_weight_wearing: float = 22.0e3
    # EC2 crack width (7.3.4) and deflection (7.4.3)
    k1: float = 1.6
    k2: float = 0.5
    k3: float = 3.4
    k4: float = 0.425
    kt: float = 0.4
    bond_xi: float = 0.6
    beta: float = 0.5
    # live load
    axle_loads: tuple = (8.0 * KIP, 32.0 * KIP, 32.0 * KIP)
    axle_spacings: tuple = (14.0 * FT, 14.0 * FT)
    impact_numerator_ft: float = 50.0
    impact_denominator_ft: float = 125.0
    impact_cap: float = 0.30
    df_divisor_ft: float = 5.5
    deflection_limit_ratio: float = 800.0


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
    "deflection_limit_ratio": "S18",
}


def field_names() -> set[str]:
    return {f.name for f in fields(Inputs)}
```

- [ ] **Step 4: Write `docs/sources.md`** (skeleton; Task 2 fills the Source/Tag/Checked columns)

```markdown
# Sources

Tags: **[primary]** opened the original; **[secondary]** an open document that quotes or reproduces the original; **[stand-in]** a chosen value, not the real CPK72 girder. Each row is checked by opening the page, not a search snippet. Quotes are kept under 15 words.

| ID | Item | Value used | Source (title, clause/table/page, link) | Tag | Checked |
|---|---|---|---|---|---|
| S01 | S01 span length | 24.50 m | candidate — verify in Task 2 | | |
| S02 | Deck layout: spacing 2.20 m, slab 0.20 m, wearing 0.05 m, 5 girders, 2 barriers at 7.0 kN/m each | as listed | candidate — verify in Task 2 | | |
| S03 | AASHTO Type IV outline (54 in deep) and A, I, yb | outline in params | candidate — verify in Task 2 | | |
| S04 | 12.7 mm seven-wire strand: area 98.7 mm², fpu 1860 MPa | as listed | candidate — verify in Task 2 | | |
| S05 | Ep for strand | 195 GPa | candidate — verify in Task 2 | | |
| S06 | Jacking 0.75 fpu; lump-sum loss 45 ksi (pretensioned) | as listed | candidate — verify in Task 2 | | |
| S07 | EC2 Table 3.1: fcm, fctm, Ecm formulas | formulas in materials.py | candidate — verify in Task 2 | | |
| S08 | EC2 §7.3.4 Eq. 7.8–7.11; k1 1.6 (tendons), k2 0.5, k3 3.4, k4 0.425, kt 0.4 | as listed | candidate — verify in Task 2 | | |
| S09 | EC2 bond ratio ξ, pretensioned strand (Table 6.2), ξ1² = ξ (Eq. 7.5) | 0.6 | candidate — verify in Task 2 | | |
| S10 | EC2 §7.4.3 Eq. 7.18–7.19; β for sustained/repeated load | 0.5 | candidate — verify in Task 2 | | |
| S11 | Unit weights: concrete 24.0, wearing surface 22.0 kN/m³ | as listed | candidate — verify in Task 2 | | |
| S12 | HS20-44 truck: 8/32/32 kip, 14 ft spacing; use by Thai DOH | as listed | candidate — verify in Task 2 | | |
| S13 | Impact I = 50/(L+125) ≤ 0.30 (L in ft) | as listed | candidate — verify in Task 2 | | |
| S14 | Live-load distribution S/5.5 wheel lines (concrete deck on PC girders) | 5.5 ft | candidate — verify in Task 2 | | |
| S15 | Strand centroid 0.10 m and cover 45 mm above soffit | as listed | candidate — verify in Task 2 | | |
| S16 | ACI 318-19 Ec = 4700√f'c, fr = 0.62√f'c (cross-check only) | formulas | candidate — verify in Task 2 | | |
| S17 | Effective flange width rule (deck on girder) | spacing governs | candidate — verify in Task 2 | | |
| S18 | Live-load deflection limit L/800 | 800 | candidate — verify in Task 2 | | |
| S19 | Concrete strengths girder 40 / slab 28 MPa | as listed | candidate — verify in Task 2 | | |
| S20 | ACI Branson Ie (318-14 §24.2.3.5) / Bischoff (318-19) | formula | candidate — verify in Task 2 | | |
| S21 | Open worked examples for tests: EC2 crack width; cracked prestressed rectangular section | fixture | candidate — verify in Task 2 | | |
| S22 | Logan §5.4 grid worked example (Plan 2) | fixture | local PDF, Plan 2 | | |
```

- [ ] **Step 5: Run to verify it passes**

Run: `python -m pytest tests/test_params.py -v`
Expected: 3 passed

- [ ] **Step 6: Commit**

```bash
git add src/crackfem/params.py docs/sources.md tests/test_params.py
git commit -m "feat: params with source ids and sources.md skeleton" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 2: B1 source pack (research, no code)

**Files:**
- Modify: `docs/sources.md`, `src/crackfem/params.py` (only if a value changes), `docs/DECISIONS.md`
- Create: `tests/fixtures/published.json`

Tools: WebSearch + WebFetch (open every page; a search snippet is not a source). Paywalled originals (EN 1992-1-1, ACI 318, Hambly) are replaced by open documents that quote the clause → tag `[secondary]`.

- [ ] **Step 1: Verify rows S01–S20.** For each row: find an open source, open it, and fill `Source` (title, clause/table/page, link), `Tag`, `Checked` (date). What to confirm per row:
  - S01: cite the as-built drawing set by sheet ("DOH BRB-CP sheet 3, approved 1995, read via the author's thesis database") — no private paths or links. Tag `[primary]`.
  - S02: tag `[stand-in]`; look for a typical concrete barrier weight per metre in an open DOT standard; if found and different from 7.0 kN/m, change `barrier_weight`.
  - S03: AASHTO Type IV dimensions (depth 54 in; flanges 20/26 in; web 8 in; haunches 6/9 in; flange thicknesses 8/8 in) **and** the published A (in²), I (in⁴), yb (in) from PCI or a state DOT girder table.
  - S04: ASTM A416 Grade 270 (1860) 0.5 in strand area 0.153 in² (98.7 mm²).
  - S05: EC2 §3.3.6(3) Ep = 195 GPa for strand.
  - S06: AASHTO Standard Specifications §9.15.1 jacking limit and Table 9.16.2.2 lump-sum loss for pretensioned members (f'c 5000 psi row).
  - S07: EC2 Table 3.1 formula column for fcm, fctm, Ecm, plus the tabulated values for C25/30, C30/37, C40/50.
  - S08: EC2 Eq. 7.8, 7.9, 7.10, 7.11 and the recommended k1, k2, k3, k4, kt.
  - S09: EC2 Table 6.2 ξ for pretensioned strand (≤ C50/60) and Eq. 7.5 "only prestressing steel → ξ1 = √ξ".
  - S10: EC2 Eq. 7.18, 7.19; β = 1.0 single short-term, 0.5 sustained or repeated.
  - S11: AASHTO LRFD Table 3.5.1-1 (or EN 1991-1-1 Annex A) unit weights.
  - S12: HS20-44 axle loads/spacing **and** an open document showing that Thai DOH bridge design uses HS20-44 (if none is found, tag the DOH-use claim `[stand-in]`).
  - S13, S14, S17, S18: AASHTO Standard Specifications §3.8.2.1, Table 3.23.1, §9.8.3 (or 8.10.1), §8.9.3.
  - S15, S19: tag `[stand-in]`; record the reasoning (Type IV bottom-flange strand rows; typical precast girder strength) and any open source found.
  - S16, S20: ACI 318-19 §19.2.2.1, §19.2.3.1; ACI 318-14 Eq. 24.2.3.5a (Branson) and the 318-19 change to Bischoff.

- [ ] **Step 2: Find open worked examples (S21).**
  - (a) One EC2 crack-width example that states σs (or Δσp), the reinforcement area and Ac,eff (or the data to compute them), c, φ, k1–k4, kt, fct,eff, Es, Ecm and the printed result wk (and, if printed, sr,max).
  - (b) One cracked **prestressed rectangular** section example (textbook or lecture notes) that states b, h, strand depth d_p, A_p, Ec, Ep, the effective prestress fpe at zero external moment, an applied moment M above cracking, and the printed neutral-axis depth and/or strand stress at M.
  Record both in `sources.md` row S21.

- [ ] **Step 3: Write `tests/fixtures/published.json`** with the values exactly as printed in the opened sources (the structure below is fixed; the numbers are what the sources print):

```json
{
  "ec2_table_3_1": {
    "source": "S07",
    "rows": [
      {"fck_MPa": 25, "fcm_MPa": 33, "fctm_MPa": 2.6, "Ecm_GPa": 31},
      {"fck_MPa": 30, "fcm_MPa": 38, "fctm_MPa": 2.9, "Ecm_GPa": 33},
      {"fck_MPa": 40, "fcm_MPa": 48, "fctm_MPa": 3.5, "Ecm_GPa": 35}
    ]
  },
  "aashto_type_iv": {"source": "S03", "h_in": 54.0, "A_in2": 789.0, "I_in4": 260741.0, "yb_in": 24.73},
  "ec2_crack_example": {
    "source": "S21",
    "inputs": {"k1": 0.8, "k2": 0.5, "k3": 3.4, "k4": 0.425, "kt": 0.4, "xi": 1.0,
               "Es": 200e9, "Ecm": 33e9, "fct_eff": 2.9e6, "cover": 0.040, "phi": 0.020,
               "As": 0.0, "Ac_eff": 0.0, "sigma_s": 0.0},
    "wk_mm": 0.0,
    "sr_max_mm": null
  },
  "pc_cracked_example": {
    "source": "S21",
    "inputs": {"b": 0.0, "h": 0.0, "d_p": 0.0, "A_p": 0.0, "Ec": 0.0, "Ep": 0.0, "fpe": 0.0, "M": 0.0},
    "neutral_axis_depth_m": null,
    "strand_stress_MPa": null
  }
}
```

The `ec2_table_3_1` and `aashto_type_iv` values shown are the candidates to confirm; replace any that the opened source prints differently. In `ec2_crack_example` and `pc_cracked_example`, every `0.0` must be replaced by the example's printed input (SI units: m, m², Pa, N·m) and at least one printed result filled in; if no open example is found, set that whole value to `null` (its test then skips with a reason) and record the gap in `sources.md` S21.

- [ ] **Step 4: Apply changes.** If any verified value differs from `params.py`, change `params.py` and add one line per change to `docs/DECISIONS.md` (`- **<field> changed from X to Y** — source Sxx, reason.`).

- [ ] **Step 5: Run the params test**

Run: `python -m pytest tests/test_params.py -v`
Expected: 3 passed

- [ ] **Step 6: Progress log entry B1** — append to `docs/PROGRESS_LOG.md`:

```markdown
## Entry 1 — B1 source pack (<date>)

**Planned:** open every source behind params.py (spec §5) and lock the values.
**Found:** <one line per row whose value changed or whose source is only [secondary]/[stand-in]>.
**Changed & why:** <params changes with source ids>.
**Numbers:** see docs/sources.md.
**Reading guide:** skim sources.md; rows tagged [stand-in] are the assumptions the examiner may question.
**Reviewed by student:** pending (after 2026-10-19).
```

(Replace each `<...>` with the actual content at the time of writing; nothing in angle brackets may remain.)

- [ ] **Step 7: Commit**

```bash
git add docs/sources.md docs/DECISIONS.md docs/PROGRESS_LOG.md tests/fixtures/published.json src/crackfem/params.py
git commit -m "docs: B1 source pack, published fixture values" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 3: `materials.py`

**Files:**
- Create: `src/crackfem/materials.py`
- Test: `tests/test_materials.py`

- [ ] **Step 1: Write the failing test**

```python
import json
from pathlib import Path

import pytest

from crackfem.materials import aci_ec, aci_fr, ec2_ecm, ec2_fcm, ec2_fctm

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("row", FIX["ec2_table_3_1"]["rows"], ids=lambda r: f"C{r['fck_MPa']}")
def test_ec2_formulas_reproduce_table_3_1(row):
    fck = row["fck_MPa"] * 1e6
    assert ec2_fcm(fck) / 1e6 == pytest.approx(row["fcm_MPa"], abs=0.01)
    assert ec2_fctm(fck) / 1e6 == pytest.approx(row["fctm_MPa"], abs=0.05)  # table rounded to 0.1 MPa
    assert ec2_ecm(fck) / 1e9 == pytest.approx(row["Ecm_GPa"], abs=0.5)     # table rounded to 1 GPa


def test_fctm_refuses_classes_above_c50():
    with pytest.raises(ValueError):
        ec2_fctm(55e6)


def test_aci_formulas_use_mpa_inside():
    assert aci_ec(40e6) == pytest.approx(4700.0 * 40.0 ** 0.5 * 1e6)
    assert aci_fr(40e6) == pytest.approx(0.62 * 40.0 ** 0.5 * 1e6)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_materials.py -v`
Expected: FAIL — `No module named 'crackfem.materials'`

- [ ] **Step 3: Write `src/crackfem/materials.py`**

```python
"""Concrete material formulas, Pa in and out.

EC2 (EN 1992-1-1 Table 3.1, source S07) is the main family. ACI 318-19
(source S16) is used only by the Branson cross-check.
"""
import math


def ec2_fcm(fck: float) -> float:
    """Mean compressive strength: fcm = fck + 8 MPa."""
    return fck + 8e6


def ec2_fctm(fck: float) -> float:
    """Mean axial tensile strength: fctm = 0.30 fck^(2/3) (MPa), valid up to C50/60."""
    if fck > 50e6:
        raise ValueError("fctm = 0.30 fck^(2/3) is valid only up to C50/60")
    return 0.30e6 * (fck / 1e6) ** (2.0 / 3.0)


def ec2_ecm(fck: float) -> float:
    """Secant modulus: Ecm = 22 (fcm / 10)^0.3 GPa, fcm in MPa."""
    return 22e9 * (ec2_fcm(fck) / 1e6 / 10.0) ** 0.3


def aci_ec(fc: float) -> float:
    """ACI 318-19 normal-weight concrete: Ec = 4700 sqrt(f'c) MPa."""
    return 4700.0 * math.sqrt(fc / 1e6) * 1e6


def aci_fr(fc: float) -> float:
    """ACI 318-19 modulus of rupture: fr = 0.62 sqrt(f'c) MPa."""
    return 0.62 * math.sqrt(fc / 1e6) * 1e6
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_materials.py -v`
Expected: all passed (3 parametrized + 2)

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/materials.py tests/test_materials.py
git commit -m "feat: EC2 and ACI concrete material formulas" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 4: `geometry.py`

**Files:**
- Create: `src/crackfem/geometry.py`
- Test: `tests/test_geometry.py`

- [ ] **Step 1: Write the failing test**

```python
import json
from pathlib import Path

import pytest

from crackfem.geometry import area_below, discretize, girder_outline, gross_properties, slab, top_of
from crackfem.params import INCH, Inputs

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))


def type_iv():
    return girder_outline(Inputs().girder_outline_in, E=1.0, fct=1.0)


def test_type_iv_matches_published_properties():
    ref = FIX["aashto_type_iv"]
    A, ybar, I = gross_properties(type_iv())
    assert top_of(type_iv()) / INCH == pytest.approx(ref["h_in"])
    assert A / INCH**2 == pytest.approx(ref["A_in2"], rel=1e-3)
    assert ybar / INCH == pytest.approx(ref["yb_in"], abs=0.02)
    assert I / INCH**4 == pytest.approx(ref["I_in4"], rel=2e-3)


def test_composite_properties_match_parallel_axis_hand_calculation():
    g = type_iv()
    Ag, yg, Ig = gross_properties(g)
    hg = top_of(g)
    n, b, t = 0.9, 2.2, 0.2
    As, ys = n * b * t, hg + t / 2
    A = Ag + As
    ybar = (Ag * yg + As * ys) / A
    I = Ig + Ag * (yg - ybar) ** 2 + n * b * t**3 / 12 + As * (ys - ybar) ** 2
    A2, y2, I2 = gross_properties(g + [slab(hg, t, b, E=n, fct=1.0)], E_ref=1.0)
    assert (A2, y2, I2) == pytest.approx((A, ybar, I), rel=1e-12)


def test_fibres_reproduce_exact_properties():
    parts = type_iv() + [slab(top_of(type_iv()), 0.2, 2.2, E=0.9, fct=1.0)]
    A, ybar, I = gross_properties(parts, E_ref=1.0)
    f = discretize(parts, dy_max=1e-3)
    A_f = (f.E * f.A).sum()
    y_f = (f.E * f.A * f.y).sum() / A_f
    I_f = (f.E * f.A * (f.y - y_f) ** 2).sum()
    assert A_f == pytest.approx(A, rel=1e-9)
    assert y_f == pytest.approx(ybar, rel=1e-7)
    assert I_f == pytest.approx(I, rel=1e-5)


def test_area_below_a_cut_inside_the_bottom_haunch():
    y_cut = 10 * INCH  # 8 in flange + 2 in into the 9 in haunch
    b_at_cut = 26 - (26 - 8) * 2 / 9
    expected_in2 = 26 * 8 + 2 * (26 + b_at_cut) / 2
    assert area_below(type_iv(), y_cut) / INCH**2 == pytest.approx(expected_in2)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_geometry.py -v`
Expected: FAIL — `No module named 'crackfem.geometry'`

- [ ] **Step 3: Write `src/crackfem/geometry.py`**

```python
"""Section outline as stacked trapezoids, exact gross properties, fibre discretisation.

Elevation y is measured upward from the girder soffit (m).
"""
import math
from dataclasses import dataclass

import numpy as np

from .params import INCH


@dataclass(frozen=True)
class Trapezoid:
    y0: float          # bottom elevation
    h: float           # height
    b0: float          # width at y0
    b1: float          # width at y0 + h
    E: float           # modulus of this concrete
    fct: float         # tensile strength of this concrete
    is_girder: bool = True

    @property
    def area(self) -> float:
        return self.h * (self.b0 + self.b1) / 2.0

    @property
    def centroid(self) -> float:
        return self.y0 + self.h * (self.b0 + 2.0 * self.b1) / (3.0 * (self.b0 + self.b1))

    @property
    def inertia_own(self) -> float:
        b0, b1 = self.b0, self.b1
        return self.h**3 * (b0 * b0 + 4.0 * b0 * b1 + b1 * b1) / (36.0 * (b0 + b1))

    def width_at(self, y):
        return self.b0 + (self.b1 - self.b0) * (y - self.y0) / self.h

    def area_below(self, y_cut: float) -> float:
        if y_cut <= self.y0:
            return 0.0
        if y_cut >= self.y0 + self.h:
            return self.area
        t = y_cut - self.y0
        return t * (self.b0 + self.width_at(y_cut)) / 2.0


def girder_outline(dims_in, E: float, fct: float) -> list[Trapezoid]:
    """Stack (height_in, width_bottom_in, width_top_in) parts from the soffit upward."""
    parts, y = [], 0.0
    for h_in, b0_in, b1_in in dims_in:
        h = h_in * INCH
        parts.append(Trapezoid(y, h, b0_in * INCH, b1_in * INCH, E, fct, True))
        y += h
    return parts


def slab(y0: float, thickness: float, width: float, E: float, fct: float) -> Trapezoid:
    return Trapezoid(y0, thickness, width, width, E, fct, False)


def top_of(parts) -> float:
    return max(p.y0 + p.h for p in parts)


def gross_properties(parts, E_ref=None) -> tuple[float, float, float]:
    """(area, centroid height, second moment about the centroid); each part is
    weighted by E / E_ref when E_ref is given (transformed section)."""
    w = [1.0 if E_ref is None else p.E / E_ref for p in parts]
    A = sum(wi * p.area for wi, p in zip(w, parts))
    ybar = sum(wi * p.area * p.centroid for wi, p in zip(w, parts)) / A
    I = sum(wi * (p.inertia_own + p.area * (p.centroid - ybar) ** 2) for wi, p in zip(w, parts))
    return A, ybar, I


def area_below(parts, y_cut: float, girder_only: bool = True) -> float:
    return sum(p.area_below(y_cut) for p in parts if p.is_girder or not girder_only)


@dataclass(frozen=True)
class Fibers:
    y: np.ndarray
    A: np.ndarray
    E: np.ndarray
    fct: np.ndarray


def discretize(parts, dy_max: float = 1e-3) -> Fibers:
    """Horizontal layers no thicker than dy_max; layer area = width at mid-height x thickness
    (exact for a linearly varying width)."""
    ys, As, Es, fs = [], [], [], []
    for p in parts:
        n = max(1, math.ceil(p.h / dy_max))
        dy = p.h / n
        y_mid = p.y0 + dy * (np.arange(n) + 0.5)
        ys.append(y_mid)
        As.append(p.width_at(y_mid) * dy)
        Es.append(np.full(n, p.E))
        fs.append(np.full(n, p.fct))
    return Fibers(np.concatenate(ys), np.concatenate(As), np.concatenate(Es), np.concatenate(fs))
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_geometry.py -v`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/geometry.py tests/test_geometry.py
git commit -m "feat: section geometry, exact properties and fibre discretisation" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 5: `section.py` — uncracked analysis, Mcr

**Files:**
- Create: `src/crackfem/section.py`
- Test: `tests/test_section_uncracked.py`

- [ ] **Step 1: Write the failing test** (independent check = classic transformed-section formulas built from `gross_properties`, not from fibres)

```python
import pytest

from crackfem.geometry import discretize, gross_properties, slab
from crackfem.section import PCSection, cracking_moment, resultants, solve_state, strand_stress, uncracked_EI

E_C, E_P = 30e9, 195e9
B, H, Y_P, A_P = 0.4, 0.8, 0.1, 1.5e-3


def rect_section(dy=2e-4):
    part = slab(0.0, H, B, E=E_C, fct=3e6)  # a rectangle is a trapezoid with b0 == b1
    return PCSection(discretize([part], dy), Y_P, A_P, E_P, H), part


def transformed():
    """Concrete rectangle + strand transformed by n = Ep/Ec (classic hand method)."""
    _, part = rect_section()
    Ac, yc, Ic = gross_properties([part])
    n = E_P / E_C
    At = Ac + n * A_P
    yt = (Ac * yc + n * A_P * Y_P) / At
    It = Ic + Ac * (yc - yt) ** 2 + n * A_P * (Y_P - yt) ** 2
    return At, yt, It


@pytest.mark.parametrize("M", [0.0, 150e3, 600e3])
def test_uncracked_state_is_in_equilibrium(M):
    sec, _ = rect_section()
    st = solve_state(sec, 4e-3, M, cracked=False)
    N, M_int = resultants(sec, 4e-3, st)
    assert abs(N) < 1e-6 * E_P * A_P * 4e-3
    assert M_int == pytest.approx(M, abs=1e-6 * 600e3)


def test_strand_stress_at_zero_moment_matches_transformed_section():
    sec, _ = rect_section()
    eps_p0 = 4e-3
    At, yt, It = transformed()
    F0 = E_P * A_P * eps_p0                     # force that the locked-in strain applies
    e = yt - Y_P
    eps_c_at_strand = -F0 / (E_C * At) - F0 * e * e / (E_C * It)
    expected = E_P * (eps_p0 + eps_c_at_strand)
    got = strand_stress(sec, eps_p0, solve_state(sec, eps_p0, 0.0, cracked=False))
    assert got == pytest.approx(expected, rel=1e-6)


def test_uncracked_EI_matches_transformed_section():
    sec, _ = rect_section()
    _, _, It = transformed()
    assert uncracked_EI(sec) == pytest.approx(E_C * It, rel=1e-6)


def test_curvature_is_linear_in_moment_with_slope_one_over_EI():
    sec, _ = rect_section()
    k0 = solve_state(sec, 4e-3, 0.0, cracked=False).kappa
    k1 = solve_state(sec, 4e-3, 300e3, cracked=False).kappa
    assert (k1 - k0) == pytest.approx(300e3 / uncracked_EI(sec), rel=1e-9)


def test_cracking_moment_matches_hand_formula():
    sec, _ = rect_section()
    eps_p0, fct = 4e-3, 3e6
    At, yt, It = transformed()
    F0 = E_P * A_P * eps_p0
    e = yt - Y_P
    # soffit stress: -F0/At - F0 e yt / It + M yt / It = fct
    Mcr = (fct + F0 / At + F0 * e * yt / It) * It / yt
    assert cracking_moment(sec, eps_p0, fct) == pytest.approx(Mcr, rel=1e-5)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_section_uncracked.py -v`
Expected: FAIL — `No module named 'crackfem.section'`

- [ ] **Step 3: Write `src/crackfem/section.py`** (cracked mode included now; Task 7 tests it)

```python
"""Elastic fibre section of a bonded prestressed member.

Conventions (shared by all of crackfem):
  * y = height above the soffit (m); tension strain/stress positive.
  * kappa > 0 is sagging; strain eps(y) = eps0 - kappa * y (eps0 = soffit strain).
  * strand stress sigma_p = Ep * (eps_p0 + eps(y_p)); eps_p0 = strand strain locked
    in when the surrounding concrete has zero strain.
  * M = external moment (sagging +). Prestress is internal: the axial resultant
    of concrete + strand is zero.
Modes: cracked=False -> concrete linear in tension and compression;
       cracked=True  -> concrete carries compression only.
"""
from dataclasses import dataclass

import numpy as np

from .geometry import Fibers


@dataclass(frozen=True)
class PCSection:
    fib: Fibers
    y_p: float
    A_p: float
    E_p: float
    h: float

    @property
    def soffit_E(self) -> float:
        return float(self.fib.E[int(np.argmin(self.fib.y))])


@dataclass(frozen=True)
class State:
    eps0: float
    kappa: float
    M: float
    cracked: bool

    def strain(self, y):
        return self.eps0 - self.kappa * y

    @property
    def neutral_axis(self) -> float:
        """Height above the soffit where the strain is zero (needs kappa != 0)."""
        return self.eps0 / self.kappa


def _solve_active(sec: PCSection, eps_p0: float, M: float, active: np.ndarray) -> tuple[float, float]:
    f = sec.fib
    ea = f.E[active] * f.A[active]
    y = f.y[active]
    epap = sec.E_p * sec.A_p
    EA = ea.sum() + epap
    ES = (ea * y).sum() + epap * sec.y_p
    EI = (ea * y * y).sum() + epap * sec.y_p**2
    J = np.array([[EA, -ES], [-ES, EI]])
    rhs = np.array([-epap * eps_p0, M + epap * sec.y_p * eps_p0])
    eps0, kappa = np.linalg.solve(J, rhs)
    return float(eps0), float(kappa)


def solve_state(sec: PCSection, eps_p0: float, M: float, cracked: bool) -> State:
    """Strain plane (eps0, kappa) that carries external moment M with zero axial force.
    Cracked mode: active-set iteration (fibres in tension drop out until the set is stable)."""
    active = np.ones(sec.fib.y.size, dtype=bool)
    eps0, kappa = _solve_active(sec, eps_p0, M, active)
    if not cracked:
        return State(eps0, kappa, M, False)
    seen = set()
    for _ in range(200):
        new_active = (eps0 - kappa * sec.fib.y) < 0.0
        if np.array_equal(new_active, active):
            return State(eps0, kappa, M, True)
        key = new_active.tobytes()
        if key in seen:
            raise RuntimeError("cracked-section active set cycled")
        seen.add(key)
        if not new_active.any():
            raise ValueError("no concrete left in compression")
        active = new_active
        eps0, kappa = _solve_active(sec, eps_p0, M, active)
    raise RuntimeError("cracked-section active set did not converge")


def concrete_stress(sec: PCSection, st: State) -> np.ndarray:
    sig = sec.fib.E * st.strain(sec.fib.y)
    return np.minimum(sig, 0.0) if st.cracked else sig


def strand_stress(sec: PCSection, eps_p0: float, st: State) -> float:
    return float(sec.E_p * (eps_p0 + st.strain(sec.y_p)))


def delta_sigma_p(sec: PCSection, st: State) -> float:
    """Strand stress change from the state of zero concrete strain at the strand level."""
    return float(sec.E_p * st.strain(sec.y_p))


def resultants(sec: PCSection, eps_p0: float, st: State) -> tuple[float, float]:
    """(axial force, internal sagging moment) of the stresses in state st."""
    sig = concrete_stress(sec, st)
    sp = strand_stress(sec, eps_p0, st)
    N = float((sig * sec.fib.A).sum() + sp * sec.A_p)
    M = float(-((sig * sec.fib.A * sec.fib.y).sum() + sp * sec.A_p * sec.y_p))
    return N, M


def uncracked_EI(sec: PCSection) -> float:
    f = sec.fib
    ea = f.E * f.A
    epap = sec.E_p * sec.A_p
    EA = ea.sum() + epap
    ES = (ea * f.y).sum() + epap * sec.y_p
    EI = (ea * f.y * f.y).sum() + epap * sec.y_p**2
    return float(EI - ES * ES / EA)


def soffit_stress(sec: PCSection, st: State) -> float:
    """Concrete stress at y = 0 for an uncracked state."""
    return sec.soffit_E * st.eps0


def cracking_moment(sec: PCSection, eps_p0: float, fct: float) -> float:
    """External moment at which the uncracked soffit stress reaches fct (fct = 0 gives the
    decompression moment). Soffit stress is linear in M, so two solves suffice."""
    s0 = soffit_stress(sec, solve_state(sec, eps_p0, 0.0, cracked=False))
    s1 = soffit_stress(sec, solve_state(sec, eps_p0, 1.0e6, cracked=False))
    return 1.0e6 * (fct - s0) / (s1 - s0)
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_section_uncracked.py -v`
Expected: 7 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/section.py tests/test_section_uncracked.py
git commit -m "feat: fibre section equilibrium, uncracked EI and cracking moment" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 6: `prestress.py`

**Files:**
- Create: `src/crackfem/prestress.py`
- Test: `tests/test_prestress.py`

- [ ] **Step 1: Write the failing test**

```python
import pytest

from crackfem.geometry import discretize, slab
from crackfem.prestress import design_prestress, effective_prestress, eps_p0_for_stress
from crackfem.section import PCSection, solve_state, soffit_stress, strand_stress

FIB = discretize([slab(0.0, 0.8, 0.4, E=30e9, fct=3e6)], 5e-4)


def make_section(A_p):
    return PCSection(FIB, 0.1, A_p, 195e9, 0.8)


def test_effective_prestress_is_jacking_minus_lump_sum_loss():
    assert effective_prestress(0.75, 1860e6, 300e6) == pytest.approx(1095e6)


def test_locked_in_strain_reproduces_the_target_stress_at_zero_moment():
    sec = make_section(1.0e-3)
    e0 = eps_p0_for_stress(sec, 1.0e9)
    assert strand_stress(sec, e0, solve_state(sec, e0, 0.0, cracked=False)) == pytest.approx(1.0e9)


def test_design_is_the_smallest_strand_count_with_no_soffit_tension():
    d = design_prestress(make_section, 98.7e-6, 1.0e9, 400e3)
    assert d.soffit_stress_service <= 0.0
    fewer = make_section((d.n_strands - 1) * 98.7e-6)
    e0 = eps_p0_for_stress(fewer, 1.0e9)
    assert soffit_stress(fewer, solve_state(fewer, e0, 400e3, cracked=False)) > 0.0
    assert d.P0 == pytest.approx(d.n_strands * 98.7e-6 * 1.0e9)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_prestress.py -v`
Expected: FAIL — `No module named 'crackfem.prestress'`

- [ ] **Step 3: Write `src/crackfem/prestress.py`**

```python
"""Prestress sizing for the stand-in girder. eta = 1 means 'as designed'.

Design rule (source S06 + spec 5.1): the smallest number of strands such that the
intact girder has no soffit tension under the full service moment.
"""
from dataclasses import dataclass
from typing import Callable

from .section import PCSection, soffit_stress, solve_state, strand_stress


@dataclass(frozen=True)
class PrestressDesign:
    n_strands: int
    A_p: float
    fpe: float                 # effective strand stress at zero external moment
    eps_p0: float              # locked-in strand strain for eta = 1
    P0: float                  # strand force at zero external moment
    M_service: float
    soffit_stress_service: float


def effective_prestress(fpj_ratio: float, fpu: float, lump_sum_loss: float) -> float:
    return fpj_ratio * fpu - lump_sum_loss


def eps_p0_for_stress(sec: PCSection, fpe: float) -> float:
    """Locked-in strain giving strand stress fpe at zero external moment (linear in eps_p0)."""
    unit = strand_stress(sec, 1.0, solve_state(sec, 1.0, 0.0, cracked=False))
    return fpe / unit


def design_prestress(make_section: Callable[[float], PCSection], strand_area: float,
                     fpe: float, M_service: float, n_max: int = 150) -> PrestressDesign:
    for n in range(1, n_max + 1):
        sec = make_section(n * strand_area)
        e0 = eps_p0_for_stress(sec, fpe)
        s = soffit_stress(sec, solve_state(sec, e0, M_service, cracked=False))
        if s <= 0.0:
            return PrestressDesign(n, n * strand_area, fpe, e0, fpe * n * strand_area, M_service, s)
    raise ValueError(f"no strand count up to {n_max} removes soffit tension")
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_prestress.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/prestress.py tests/test_prestress.py
git commit -m "feat: zero-tension prestress sizing" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 7: `section.py` — cracked analysis tests

**Files:**
- Test: `tests/test_section_cracked.py` (no production change expected; fix `section.py` only if a test fails)

- [ ] **Step 1: Write the test** (independent check = textbook cracked RC section: kd = d(√(2ρn + (ρn)²) − ρn), Icr = b(kd)³/3 + nAs(d − kd)²)

```python
import json
import math
from pathlib import Path

import numpy as np
import pytest

from crackfem.geometry import Trapezoid, discretize, slab, top_of
from crackfem.section import (PCSection, concrete_stress, delta_sigma_p, resultants, solve_state,
                              strand_stress)

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))
E_C, E_S = 25e9, 200e9
B, H, D, A_S = 0.3, 0.6, 0.55, 1.5e-3


def rc_section():
    return PCSection(discretize([slab(0.0, H, B, E=E_C, fct=3e6)], 1e-4), H - D, A_S, E_S, H)


def test_cracked_rc_section_matches_textbook_formulas():
    n = E_S / E_C
    rho = A_S / (B * D)
    kd = D * (math.sqrt(2 * rho * n + (rho * n) ** 2) - rho * n)
    Icr = B * kd**3 / 3 + n * A_S * (D - kd) ** 2
    M = 100e3
    st = solve_state(rc_section(), 0.0, M, cracked=True)
    assert H - st.neutral_axis == pytest.approx(kd, abs=2e-4)
    assert st.kappa == pytest.approx(M / (E_C * Icr), rel=2e-3)


def pc_section():
    return PCSection(discretize([slab(0.0, 0.8, 0.4, E=30e9, fct=3e6)], 5e-4), 0.1, 1.5e-3, 195e9, 0.8)


@pytest.mark.parametrize("M", [100e3, 400e3, 700e3])
def test_cracked_pc_state_is_in_equilibrium_with_no_concrete_tension(M):
    sec = pc_section()
    st = solve_state(sec, 4e-3, M, cracked=True)
    N, M_int = resultants(sec, 4e-3, st)
    assert abs(N) < 1e-6 * 195e9 * 1.5e-3 * 4e-3
    assert M_int == pytest.approx(M, rel=1e-6)
    assert concrete_stress(sec, st).max() <= 0.0


def test_below_decompression_the_cracked_state_equals_the_uncracked_state():
    sec = pc_section()
    a = solve_state(sec, 4e-3, 50e3, cracked=False)
    b = solve_state(sec, 4e-3, 50e3, cracked=True)
    assert (b.eps0, b.kappa) == pytest.approx((a.eps0, a.kappa), rel=1e-12)


def test_strand_stress_change_is_positive_only_when_the_neutral_axis_is_above_the_strand():
    sec = pc_section()
    for M in np.linspace(50e3, 800e3, 16):
        st = solve_state(sec, 4e-3, M, cracked=True)
        ds = delta_sigma_p(sec, st)
        if ds > 0:
            assert st.kappa > 0 and st.neutral_axis > sec.y_p


def test_published_cracked_prestressed_example():
    ex = FIX.get("pc_cracked_example")
    if ex is None:
        pytest.skip("no open cracked prestressed-section example located (docs/sources.md S21)")
    i = ex["inputs"]
    parts = [Trapezoid(y0, h, b0, b1, i["Ec"], 0.0) for y0, h, b0, b1 in i["parts"]]
    h_tot = top_of(parts)
    sec = PCSection(discretize(parts, 1e-4), h_tot - i["d_p"], i["A_p"], i["Ep"], h_tot)
    e0 = i["fdc"] / i["Ep"]   # locked-in strain = strand strain at zero concrete strain = fdc / Ep
    st = solve_state(sec, e0, i["M"], cracked=True)
    if ex["neutral_axis_depth_m"] is not None:
        assert h_tot - st.neutral_axis == pytest.approx(ex["neutral_axis_depth_m"], rel=0.03)
    if ex["strand_stress_MPa"] is not None:
        assert strand_stress(sec, e0, st) / 1e6 == pytest.approx(ex["strand_stress_MPa"], rel=0.03)
```

- [ ] **Step 2: Run**

Run: `python -m pytest tests/test_section_cracked.py -v`
Expected: 7 passed, or 6 passed + 1 skipped (S21 example not found). If `test_cracked_rc_section_matches_textbook_formulas` fails, the fault is in `solve_state` cracked mode (check the active-set test `< 0.0` and the sign of `rhs`) — fix `section.py`, not the test.

- [ ] **Step 3: Commit**

```bash
git add tests/test_section_cracked.py src/crackfem/section.py
git commit -m "test: cracked section against textbook RC formulas and equilibrium" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 8: `crackwidth.py`

**Files:**
- Create: `src/crackfem/crackwidth.py`
- Test: `tests/test_crackwidth.py`

- [ ] **Step 1: Write the failing test**

```python
import json
from pathlib import Path

import pytest

from crackfem.crackwidth import EC2Crack, effective_tension_height

FIX = json.loads((Path(__file__).parent / "fixtures" / "published.json").read_text(encoding="utf-8"))

PC = EC2Crack(k1=1.6, k2=0.5, k3=3.4, k4=0.425, kt=0.4, xi=0.6, Ep=195e9, Ecm=35e9,
              fct_eff=3.5e6, cover=0.045, phi=12.7e-3, A_p=2.6e-3)
AC = 0.13


def test_published_ec2_example():
    ex = FIX["ec2_crack_example"]
    if ex is None:
        pytest.skip("no open EC2 worked example located (docs/sources.md S21)")
    i = ex["inputs"]
    c = EC2Crack(i["k1"], i["k2"], i["k3"], i["k4"], i["kt"], i["xi"], i["Es"], i["Ecm"],
                 i["fct_eff"], i["cover"], i["phi"], i["As"])
    assert c.width(i["sigma_s"], i["Ac_eff"]) * 1e3 == pytest.approx(ex["wk_mm"], abs=0.005)
    if ex.get("sr_max_mm") is not None:
        assert c.sr_max(i["Ac_eff"]) * 1e3 == pytest.approx(ex["sr_max_mm"], rel=0.01)


def test_no_width_without_tension_in_the_strand():
    assert PC.width(0.0, AC) == 0.0
    assert PC.width(-50e6, AC) == 0.0


@pytest.mark.parametrize("ds", [10e6, 60e6, 200e6, 450e6])
def test_width_to_stress_round_trip_on_both_branches(ds):
    assert PC.stress_from_width(PC.width(ds, AC), AC) == pytest.approx(ds, rel=1e-12)


def test_branches_meet_at_two_and_a_half_times_the_offset():
    rho = PC.rho(AC)
    C = PC.kt * PC.fct_eff / rho * (1 + PC.Ep / PC.Ecm * rho)
    ds = 2.5 * C
    assert (ds - C) / PC.Ep == pytest.approx(0.6 * ds / PC.Ep, rel=1e-12)
    assert PC.width(ds * (1 - 1e-9), AC) == pytest.approx(PC.width(ds * (1 + 1e-9), AC), rel=1e-6)


def test_effective_tension_height_takes_the_smallest_limit():
    assert effective_tension_height(h=1.5, y_p=0.1, y_na=0.9) == pytest.approx(0.25)   # 2.5 y_p
    assert effective_tension_height(h=1.5, y_p=0.2, y_na=0.6) == pytest.approx(0.20)   # y_na / 3
    assert effective_tension_height(h=0.4, y_p=0.1, y_na=0.9) == pytest.approx(0.20)   # h / 2
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_crackwidth.py -v`
Expected: FAIL — `No module named 'crackfem.crackwidth'`

- [ ] **Step 3: Write `src/crackfem/crackwidth.py`**

```python
"""EC2 (EN 1992-1-1) section 7.3.4 crack width for bonded prestressed members (S08, S09).

  wk = sr,max * (eps_sm - eps_cm)                                           (7.8)
  eps_sm - eps_cm = [ds - kt fct,eff / rho (1 + alpha_e rho)] / Ep >= 0.6 ds / Ep  (7.9)
  rho = rho_p,eff = xi1^2 A_p / Ac,eff, with xi1^2 = xi for prestressing steel only (7.10, 7.5)
  sr,max = k3 c + k1 k2 k4 phi / rho                                         (7.11)
ds = strand stress change from the state of zero concrete strain at the strand level
(delta sigma_p); it replaces sigma_s for prestressed members. For ordinary bars set
xi = 1 and pass the bar area as A_p (used by the published-example test).
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class EC2Crack:
    k1: float
    k2: float
    k3: float
    k4: float
    kt: float
    xi: float
    Ep: float
    Ecm: float
    fct_eff: float
    cover: float
    phi: float
    A_p: float

    def rho(self, Ac_eff: float) -> float:
        return self.xi * self.A_p / Ac_eff

    def sr_max(self, Ac_eff: float) -> float:
        return self.k3 * self.cover + self.k1 * self.k2 * self.k4 * self.phi / self.rho(Ac_eff)

    def _offset(self, Ac_eff: float) -> float:
        r = self.rho(Ac_eff)
        return self.kt * self.fct_eff / r * (1.0 + self.Ep / self.Ecm * r)

    def strain_difference(self, ds: float, Ac_eff: float) -> float:
        if ds <= 0.0:
            return 0.0
        return max((ds - self._offset(Ac_eff)) / self.Ep, 0.6 * ds / self.Ep)

    def width(self, ds: float, Ac_eff: float) -> float:
        return self.sr_max(Ac_eff) * self.strain_difference(ds, Ac_eff)

    def stress_from_width(self, w: float, Ac_eff: float) -> float:
        """Closed-form inverse of width(); the 0.6 floor governs while ds <= 2.5 * offset."""
        if w <= 0.0:
            return 0.0
        sr = self.sr_max(Ac_eff)
        C = self._offset(Ac_eff)
        if w <= sr * 1.5 * C / self.Ep:
            return w * self.Ep / (0.6 * sr)
        return w * self.Ep / sr + C


def effective_tension_height(h: float, y_p: float, y_na: float) -> float:
    """hc,ef = min(2.5 (h - d), (h - x) / 3, h / 2), with h - d = y_p and h - x = y_na."""
    return min(2.5 * y_p, y_na / 3.0, h / 2.0)
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_crackwidth.py -v`
Expected: all passed (the published-example test passes or skips with the S21 reason)

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/crackwidth.py tests/test_crackwidth.py
git commit -m "feat: EC2 crack width for prestressed members with closed-form inverse" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 9: `stiffness.py`

**Files:**
- Create: `src/crackfem/stiffness.py`
- Test: `tests/test_stiffness.py`

Staged construction (decision 2026-10-12): the FE model carries only the loads applied after a reference state with total moment `M_ref` (`M_ref = 0` for the prestressed, unloaded state; `M_ref` = the stage-A moment for a staged girder). The secant stiffness is therefore `EI = (M − M_ref) / (κ(M) − κ(M_ref))`, with `M` and `M_ref` total moments and κ the EC2 7.4.3 curvature.

- [ ] **Step 1: Write the failing test**

```python
import numpy as np
import pytest

from crackfem.geometry import discretize, slab
from crackfem.section import PCSection, cracking_moment, solve_state, uncracked_EI
from crackfem.stiffness import branson_secant_EI, ec2_secant_EI

SEC = PCSection(discretize([slab(0.0, 0.8, 0.4, E=30e9, fct=3e6)], 5e-4), 0.1, 1.5e-3, 195e9, 0.8)
E0 = 4e-3
MCR = cracking_moment(SEC, E0, 3e6)


def kI(M):
    return solve_state(SEC, E0, M, False).kappa


def kII(M):
    return solve_state(SEC, E0, M, True).kappa


def test_below_cracking_the_secant_stiffness_is_the_uncracked_stiffness():
    for M_ref in (0.0, 0.5 * MCR):
        assert ec2_secant_EI(SEC, E0, 0.9 * MCR, MCR, 0.5, M_ref) == pytest.approx(uncracked_EI(SEC), rel=1e-9)
        assert branson_secant_EI(SEC, E0, 0.9 * MCR, MCR, M_ref) == pytest.approx(uncracked_EI(SEC), rel=1e-9)


def test_ec2_interpolation_with_beta_one_at_twice_mcr():
    M = 2 * MCR
    expected = M / (0.75 * kII(M) + 0.25 * kI(M) - kI(0.0))      # zeta = 1 - (1/2)^2
    assert ec2_secant_EI(SEC, E0, M, MCR, 1.0) == pytest.approx(expected, rel=1e-12)


def test_secant_from_a_reference_moment_is_increment_over_increment():
    M, M_ref = 2 * MCR, 0.5 * MCR
    expected = (M - M_ref) / (0.75 * kII(M) + 0.25 * kI(M) - kI(M_ref))
    assert ec2_secant_EI(SEC, E0, M, MCR, 1.0, M_ref) == pytest.approx(expected, rel=1e-12)


def test_no_load_increment_keeps_the_uncracked_stiffness():
    assert ec2_secant_EI(SEC, E0, 1.5 * MCR, MCR, 0.5, M_ref=1.5 * MCR) == pytest.approx(uncracked_EI(SEC))


def test_secant_stiffness_falls_as_moment_grows_above_cracking():
    Ms = np.linspace(1.05 * MCR, 2.5 * MCR, 8)
    for M_ref in (0.0, 0.5 * MCR):
        ei = [ec2_secant_EI(SEC, E0, m, MCR, 0.5, M_ref) for m in Ms]
        assert all(b < a for a, b in zip(ei, ei[1:]))
        assert ei[0] < uncracked_EI(SEC)


def test_branson_lies_between_cracked_and_uncracked():
    M = 2 * MCR
    EI_cr = M / (kII(M) - kI(0.0))
    assert EI_cr < branson_secant_EI(SEC, E0, M, MCR) < uncracked_EI(SEC)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_stiffness.py -v`
Expected: FAIL — `No module named 'crackfem.stiffness'`

- [ ] **Step 3: Write `src/crackfem/stiffness.py`**

```python
"""Secant flexural stiffness of one element from its moment (S10 main, S20 cross-check).

The FE model carries only the loads applied after a reference state with total moment M_ref
(0 = prestressed and unloaded; the stage-A moment for a staged girder), so
EI_sec = (M - M_ref) / (kappa(M) - kappa(M_ref)), with M and M_ref total moments.
Hogging (M <= 0) keeps the uncracked stiffness: the girders are simply supported and top
cracking is not modelled.
"""
from .section import PCSection, solve_state, uncracked_EI


def _kappa_uncracked(sec: PCSection, eps_p0: float, M: float) -> float:
    return solve_state(sec, eps_p0, M, cracked=False).kappa


def ec2_kappa(sec: PCSection, eps_p0: float, M: float, Mcr: float, beta: float) -> float:
    """EC2 7.4.3: kappa = zeta kappa_II + (1 - zeta) kappa_I, zeta = 1 - beta (Mcr / M)^2."""
    Mcr = max(Mcr, 0.0)
    kI = _kappa_uncracked(sec, eps_p0, M)
    if M <= 0.0 or M <= Mcr:
        return kI
    kII = solve_state(sec, eps_p0, M, cracked=True).kappa
    zeta = 1.0 - beta * (Mcr / M) ** 2
    return zeta * kII + (1.0 - zeta) * kI


def branson_kappa(sec: PCSection, eps_p0: float, M: float, Mcr: float) -> float:
    """ACI Branson form on secant stiffness from the zero-moment state:
    kappa = kappa(0) + M / (r EI_I + (1 - r) EI_cr), r = (Mcr / M)^3."""
    Mcr = max(Mcr, 0.0)
    if M <= 0.0 or M <= Mcr:
        return _kappa_uncracked(sec, eps_p0, M)
    k0 = _kappa_uncracked(sec, eps_p0, 0.0)
    EI_cr = M / (solve_state(sec, eps_p0, M, cracked=True).kappa - k0)
    r = (Mcr / M) ** 3
    return k0 + M / (r * uncracked_EI(sec) + (1.0 - r) * EI_cr)


def ec2_secant_EI(sec: PCSection, eps_p0: float, M: float, Mcr: float, beta: float,
                  M_ref: float = 0.0) -> float:
    if M <= M_ref or M <= max(Mcr, 0.0):
        return uncracked_EI(sec)
    return (M - M_ref) / (ec2_kappa(sec, eps_p0, M, Mcr, beta) - ec2_kappa(sec, eps_p0, M_ref, Mcr, beta))


def branson_secant_EI(sec: PCSection, eps_p0: float, M: float, Mcr: float, M_ref: float = 0.0) -> float:
    if M <= M_ref or M <= max(Mcr, 0.0):
        return uncracked_EI(sec)
    return (M - M_ref) / (branson_kappa(sec, eps_p0, M, Mcr) - branson_kappa(sec, eps_p0, M_ref, Mcr))
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_stiffness.py -v`
Expected: 6 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/stiffness.py tests/test_stiffness.py
git commit -m "feat: EC2 and Branson secant stiffness from a reference moment" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---


## Task 10: `elements.py` (beam element, symbolic check)

**Files:**
- Create: `src/crackfem/elements.py`
- Test: `tests/test_elements.py`

- [ ] **Step 1: Write the failing test** (independent evidence = sympy energy/work integrals of the Hermite shape functions)

```python
import numpy as np
import sympy as sp

from crackfem.elements import beam_k, hermite, point_equivalent, udl_equivalent

x, L, EI, q, P, a = sp.symbols("x L EI q P a", positive=True)
xi = x / L
N = sp.Matrix([1 - 3 * xi**2 + 2 * xi**3,
               L * (xi - 2 * xi**2 + xi**3),
               3 * xi**2 - 2 * xi**3,
               L * (-xi**2 + xi**3)])
CASES = [(2.0e9, 1.3), (5.5, 0.4)]


def test_shape_functions_interpolate_nodal_values_and_slopes():
    dN = N.diff(x)
    assert list(N.subs(x, 0)) == [1, 0, 0, 0]
    assert list(dN.subs(x, 0)) == [0, 1, 0, 0]
    assert list(N.subs(x, L)) == [0, 0, 1, 0]
    assert list(sp.simplify(dN.subs(x, L))) == [0, 0, 0, 1]


def test_stiffness_equals_bending_energy_integral():
    B = N.diff(x, 2)
    k_sym = sp.lambdify((EI, L), sp.integrate(EI * B * B.T, (x, 0, L)))
    for EIv, Lv in CASES:
        np.testing.assert_allclose(beam_k(EIv, Lv), np.array(k_sym(EIv, Lv), float), rtol=1e-12)


def test_udl_equivalent_equals_work_integral():
    f_sym = sp.lambdify((q, L), sp.integrate(q * N, (x, 0, L)))
    np.testing.assert_allclose(udl_equivalent(-12.5, 1.7), np.array(f_sym(-12.5, 1.7), float).ravel(), rtol=1e-12)


def test_point_equivalent_equals_shape_functions_at_the_load():
    f_sym = sp.lambdify((P, a, L), P * N.subs(x, a))
    np.testing.assert_allclose(point_equivalent(-3.0, 0.45, 1.2),
                               np.array(f_sym(-3.0, 0.45, 1.2), float).ravel(), rtol=1e-12)
    np.testing.assert_allclose(hermite(0.45 / 1.2, 1.2),
                               np.array(sp.lambdify((a, L), N.subs(x, a))(0.45, 1.2), float).ravel(), rtol=1e-12)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_elements.py -v`
Expected: FAIL — `No module named 'crackfem.elements'`

- [ ] **Step 3: Write `src/crackfem/elements.py`**

```python
"""Euler-Bernoulli beam element (course Ch3). DOFs per element: (v1, theta1, v2, theta2);
v up positive, theta counter-clockwise positive, loads up positive."""
import numpy as np


def beam_k(EI: float, L: float) -> np.ndarray:
    c = EI / L**3
    return c * np.array([
        [12.0, 6.0 * L, -12.0, 6.0 * L],
        [6.0 * L, 4.0 * L * L, -6.0 * L, 2.0 * L * L],
        [-12.0, -6.0 * L, 12.0, -6.0 * L],
        [6.0 * L, 2.0 * L * L, -6.0 * L, 4.0 * L * L],
    ])


def hermite(xi: float, L: float) -> np.ndarray:
    """Cubic Hermite shape functions at xi = x / L."""
    return np.array([1 - 3 * xi**2 + 2 * xi**3,
                     L * (xi - 2 * xi**2 + xi**3),
                     3 * xi**2 - 2 * xi**3,
                     L * (-xi**2 + xi**3)])


def udl_equivalent(q: float, L: float) -> np.ndarray:
    """Consistent nodal loads of a uniform load q (up +) over the element."""
    return np.array([q * L / 2.0, q * L * L / 12.0, q * L / 2.0, -q * L * L / 12.0])


def point_equivalent(P: float, a: float, L: float) -> np.ndarray:
    """Consistent nodal loads of a point load P (up +) at distance a from node 1."""
    return P * hermite(a / L, L)
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_elements.py -v`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/elements.py tests/test_elements.py
git commit -m "feat: beam element with symbolically verified stiffness and loads" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 11: `model.py` + `solve.py`

**Files:**
- Create: `src/crackfem/model.py`, `src/crackfem/solve.py`
- Test: `tests/test_solve.py`

- [ ] **Step 1: Write the failing test**

```python
import numpy as np
import pytest
from scipy.integrate import quad

from crackfem.model import BeamModel, simply_supported, with_loads
from crackfem.solve import midpoint_moments, moment_at, node_moments, solve

L, EI, W = 24.5, 7.5e9, 28e3


@pytest.mark.parametrize("n", [1, 2, 3, 5, 8, 40])
def test_simply_supported_udl_is_exact_at_the_nodes_for_any_mesh(n):
    m = with_loads(simply_supported(L, n), q=-W)
    sol = solve(m, EI)
    x = m.x
    exact = -W * x * (L**3 - 2 * L * x**2 + x**3) / (24 * EI)
    np.testing.assert_allclose(sol.D[0::2], exact, rtol=1e-10, atol=1e-15)


def test_reactions_balance_the_applied_loads():
    m = with_loads(simply_supported(L, 10), q=-W, point_loads=[(3.3, -100e3), (12.25, -50e3)])
    sol = solve(m, EI)
    assert sol.R.sum() == pytest.approx(W * L + 150e3, rel=1e-12)


def test_cantilever_tip_load():
    m = BeamModel(x=np.linspace(0, 3.0, 7), fixed_dofs=(0, 1))
    m = with_loads(m, q=0.0, point_loads=[(3.0, 10e3)])
    sol = solve(m, 2.0e7)
    assert sol.D[-2] == pytest.approx(10e3 * 3.0**3 / (3 * 2.0e7), rel=1e-12)


def test_propped_cantilever_fixed_end_moment():
    m = BeamModel(x=np.linspace(0, L, 9), fixed_dofs=(0, 1, 16))
    m = with_loads(m, q=-W)
    sol = solve(m, EI)
    assert moment_at(m, sol, 0, 0.0) == pytest.approx(-W * L**2 / 8, rel=1e-10)


def test_stepped_stiffness_matches_virtual_work():
    n = 10
    m = with_loads(simply_supported(L, n), q=-W)
    EI_el = np.where(m.midpoints < L / 2, EI, 0.6 * EI)
    sol = solve(m, EI_el)
    M = lambda x: W * x * (L - x) / 2
    unit = lambda x: x / 2 if x <= L / 2 else (L - x) / 2
    d1 = quad(lambda x: M(x) * unit(x) / EI, 0, L / 2)[0]
    d2 = quad(lambda x: M(x) * unit(x) / (0.6 * EI), L / 2, L)[0]
    assert -sol.D[2 * (n // 2)] == pytest.approx(d1 + d2, rel=1e-9)


def test_midpoint_moments_match_statics_with_point_loads():
    pts = [(5.0, -80e3), (17.0, -40e3)]
    m = with_loads(simply_supported(L, 12), q=-W, point_loads=pts)
    sol = solve(m, EI)
    RA = W * L / 2 + sum(-P * (L - xp) / L for xp, P in pts)
    def statics(x):
        return RA * x - W * x * x / 2 - sum(-P * (x - xp) for xp, P in pts if xp < x)
    np.testing.assert_allclose(midpoint_moments(m, sol), [statics(x) for x in m.midpoints], rtol=1e-9)


def test_node_moments_match_statics():
    pts = [(5.0, -80e3), (17.0, -40e3)]
    m = with_loads(simply_supported(L, 12), q=-W, point_loads=pts)
    sol = solve(m, EI)
    RA = W * L / 2 + sum(-P * (L - xp) / L for xp, P in pts)
    statics = [RA * x - W * x * x / 2 - sum(-P * (x - xp) for xp, P in pts if xp < x) for x in m.x]
    np.testing.assert_allclose(node_moments(m, sol), statics, rtol=1e-9, atol=1e-6)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_solve.py -v`
Expected: FAIL — `No module named 'crackfem.model'`

- [ ] **Step 3: Write `src/crackfem/model.py`**

```python
"""1-D beam model data. DOFs per node (v, theta) -> global (2i, 2i+1). Units N, m."""
from dataclasses import dataclass, field, replace

import numpy as np


@dataclass(frozen=True)
class BeamModel:
    x: np.ndarray                                   # node coordinates
    fixed_dofs: tuple                               # restrained global DOFs
    q: np.ndarray | None = field(default=None)      # per-element uniform load, up +
    point_loads: tuple = ()                         # ((x, P), ...), up +

    @property
    def n_el(self) -> int:
        return self.x.size - 1

    @property
    def n_dof(self) -> int:
        return 2 * self.x.size

    @property
    def lengths(self) -> np.ndarray:
        return np.diff(self.x)

    @property
    def midpoints(self) -> np.ndarray:
        return 0.5 * (self.x[:-1] + self.x[1:])

    def udl(self) -> np.ndarray:
        return np.zeros(self.n_el) if self.q is None else self.q

    def locate(self, xp: float) -> tuple[int, float]:
        """(element index, distance from its first node) of a point on the beam."""
        if xp < self.x[0] or xp > self.x[-1]:
            raise ValueError(f"x = {xp} is outside the beam")
        e = min(int(np.searchsorted(self.x, xp, side="right")) - 1, self.n_el - 1)
        return e, float(xp - self.x[e])


def simply_supported(L: float, n_el: int) -> BeamModel:
    return BeamModel(x=np.linspace(0.0, L, n_el + 1), fixed_dofs=(0, 2 * n_el))


def with_loads(model: BeamModel, q: float = 0.0, point_loads=()) -> BeamModel:
    return replace(model, q=np.full(model.n_el, float(q)), point_loads=tuple(point_loads))
```

- [ ] **Step 4: Write `src/crackfem/solve.py`**

```python
"""Direct stiffness assembly and partition solve (course method)."""
from dataclasses import dataclass

import numpy as np

from .elements import beam_k, point_equivalent, udl_equivalent
from .model import BeamModel


@dataclass(frozen=True)
class Solution:
    D: np.ndarray             # nodal displacements (v, theta per node)
    R: np.ndarray             # reactions at restrained DOFs, zero elsewhere
    end_forces: np.ndarray    # (n_el, 4): [V1, M1, V2, M2] acting on each element


def _dofs(e: int) -> np.ndarray:
    return np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3])


def element_loads(model: BeamModel) -> list[np.ndarray]:
    q, L = model.udl(), model.lengths
    f = [udl_equivalent(q[e], L[e]) for e in range(model.n_el)]
    for xp, P in model.point_loads:
        e, a = model.locate(xp)
        f[e] = f[e] + point_equivalent(P, a, L[e])
    return f


def assemble(model: BeamModel, EI: np.ndarray):
    K = np.zeros((model.n_dof, model.n_dof))
    F = np.zeros(model.n_dof)
    fl = element_loads(model)
    for e in range(model.n_el):
        idx = _dofs(e)
        K[np.ix_(idx, idx)] += beam_k(EI[e], model.lengths[e])
        F[idx] += fl[e]
    return K, F, fl


def solve(model: BeamModel, EI) -> Solution:
    EI = np.broadcast_to(np.asarray(EI, dtype=float), (model.n_el,))
    K, F, fl = assemble(model, EI)
    fixed = np.array(model.fixed_dofs)
    free = np.setdiff1d(np.arange(model.n_dof), fixed)
    D = np.zeros(model.n_dof)
    D[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])   # partition: K_ff D_f = F_f
    R = np.zeros(model.n_dof)
    R[fixed] = K[fixed] @ D - F[fixed]
    ef = np.array([beam_k(EI[e], model.lengths[e]) @ D[_dofs(e)] - fl[e] for e in range(model.n_el)])
    return Solution(D, R, ef)


def moment_at(model: BeamModel, sol: Solution, e: int, xl: float) -> float:
    """Internal sagging moment at distance xl from the first node of element e."""
    V1, M1 = sol.end_forces[e, 0], sol.end_forces[e, 1]
    m = -M1 + V1 * xl + model.udl()[e] * xl * xl / 2.0
    for xp, P in model.point_loads:
        ee, a = model.locate(xp)
        if ee == e and a < xl:
            m += P * (xl - a)
    return float(m)


def midpoint_moments(model: BeamModel, sol: Solution) -> np.ndarray:
    return np.array([moment_at(model, sol, e, model.lengths[e] / 2.0) for e in range(model.n_el)])


def node_moments(model: BeamModel, sol: Solution) -> np.ndarray:
    """Internal sagging moment at every node (start of each element, plus the last node)."""
    m = [moment_at(model, sol, e, 0.0) for e in range(model.n_el)]
    m.append(moment_at(model, sol, model.n_el - 1, model.lengths[-1]))
    return np.array(m)
```

- [ ] **Step 5: Run to verify it passes**

Run: `python -m pytest tests/test_solve.py -v`
Expected: 12 passed (6 parametrized + 6)

- [ ] **Step 6: Commit**

```bash
git add src/crackfem/model.py src/crackfem/solve.py tests/test_solve.py
git commit -m "feat: beam model, assembly, partition solve and internal moments" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 12: `loads.py`

**Files:**
- Create: `src/crackfem/loads.py`
- Test: `tests/test_loads.py`

- [ ] **Step 0: Add the DOH live-load factor parameter** (controller decision 2026-10-07, source S12: JICA 2000 Table 3.1 "Main girder : DOH standard (AASHTO HS 20-44x1.30)")

In `src/crackfem/params.py`, add to `Inputs` right after `df_divisor_ft`:

```python
    live_load_factor: float = 1.30            # DOH design truck = HS20-44 x 1.30 (S12)
```

and to `SOURCES`: `"live_load_factor": "S12",`. In `docs/sources.md` row S12, append to the "Value used" cell: `; live-load factor 1.30 (DOH standard HS20-44 x 1.30)`. Run `python -m pytest tests/test_params.py -q` → 3 passed.

- [ ] **Step 1: Write the failing test** (independent check = for a triangular influence line the maximum is reached with an axle at the apex, so only 3 placements × 2 directions need hand evaluation)

```python
import pytest

from crackfem.loads import (axle_offsets, dead_load_line, impact_factor, midspan_moment_statics,
                            service_truck_loads, stage_a_line, stage_b_line, truck_fraction,
                            truck_position_max_midspan)
from crackfem.params import FT, KIP, Inputs

INP = Inputs()


def test_axle_loads_in_newtons():
    assert [P for _, P in axle_offsets(INP)] == pytest.approx([8 * KIP, 32 * KIP, 32 * KIP])
    assert axle_offsets(INP, reverse=True)[0][1] == pytest.approx(32 * KIP)


def test_impact_factor_for_s01():
    assert impact_factor(INP) == pytest.approx(50 / (24.5 / FT + 125))


def test_truck_fraction_is_half_of_s_over_5_5():
    assert truck_fraction(INP) == pytest.approx((2.2 / FT) / 5.5 / 2)


def test_scan_finds_the_apex_placement_maximum():
    L = INP.span
    il = lambda x: x / 2 if x <= L / 2 else (L - x) / 2
    best = 0.0
    for rev in (False, True):
        axles = axle_offsets(INP, rev)
        for o_k, _ in axles:
            x0 = L / 2 - o_k
            best = max(best, sum(P * il(x0 + o) for o, P in axles if 0 <= x0 + o <= L))
    _, m = truck_position_max_midspan(INP)
    assert m == pytest.approx(best, rel=1e-6)


def test_dead_load_is_stage_a_plus_stage_b():
    A = 0.5
    a = INP.unit_weight_concrete * (A + 0.2 * 2.2)
    b = INP.unit_weight_wearing * 0.05 * 2.2 + 2 * INP.barrier_weight / 5
    assert stage_a_line(INP, A) == pytest.approx(a)
    assert stage_b_line(INP) == pytest.approx(b)
    assert dead_load_line(INP, A) == pytest.approx(a + b)


def test_service_moment_combines_dead_and_scaled_truck():
    pos, m_truck = truck_position_max_midspan(INP)
    scale = truck_fraction(INP) * (1 + impact_factor(INP)) * INP.live_load_factor
    M = midspan_moment_statics(INP.span, 20e3, service_truck_loads(INP))
    assert M == pytest.approx(20e3 * INP.span**2 / 8 + scale * m_truck, rel=1e-9)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_loads.py -v`
Expected: FAIL — `No module named 'crackfem.loads'`

- [ ] **Step 3: Write `src/crackfem/loads.py`**

```python
"""Dead load per girder and the HS20-44 design truck (S02, S11, S12, S13, S14).
Magnitudes are positive here; the caller applies them downward."""
import numpy as np

from .params import FT, Inputs


def stage_a_line(inp: Inputs, girder_area: float) -> float:
    """Girder + wet slab per girder, carried by the bare girder (stage A)."""
    return inp.unit_weight_concrete * (girder_area + inp.slab_thickness * inp.girder_spacing)


def stage_b_line(inp: Inputs) -> float:
    """Wearing surface + barrier share per girder, carried by the composite girder (stage B)."""
    return (inp.unit_weight_wearing * inp.wearing_thickness * inp.girder_spacing
            + inp.n_barriers * inp.barrier_weight / inp.n_girders)


def dead_load_line(inp: Inputs, girder_area: float) -> float:
    return stage_a_line(inp, girder_area) + stage_b_line(inp)


def impact_factor(inp: Inputs) -> float:
    return min(inp.impact_numerator_ft / (inp.span / FT + inp.impact_denominator_ft), inp.impact_cap)


def truck_fraction(inp: Inputs) -> float:
    """Fraction of one truck per girder: (S / 5.5) wheel lines divided by 2 wheel lines."""
    return (inp.girder_spacing / FT) / inp.df_divisor_ft / 2.0


def axle_offsets(inp: Inputs, reverse: bool = False) -> list[tuple[float, float]]:
    offs = np.concatenate([[0.0], np.cumsum(inp.axle_spacings)])
    loads = list(inp.axle_loads)
    if reverse:
        offs, loads = offs[-1] - offs[::-1], loads[::-1]
    return list(zip(offs.tolist(), loads))


def truck_position_max_midspan(inp: Inputs, step: float = 1e-3):
    """Axle positions on the span that maximise the midspan moment, and that moment."""
    L = inp.span
    best_pos, best_m = None, -1.0
    for rev in (False, True):
        axles = axle_offsets(inp, rev)
        offs = np.array([o for o, _ in axles])
        loads = np.array([P for _, P in axles])
        x0 = np.arange(-offs[-1], L + step / 2, step)
        X = x0[:, None] + offs[None, :]
        on = (X >= 0.0) & (X <= L)
        il = np.where(X <= L / 2, X / 2, (L - X) / 2)
        m = (np.where(on, il, 0.0) * loads[None, :]).sum(axis=1)
        i = int(np.argmax(m))
        if m[i] > best_m:
            best_m = float(m[i])
            best_pos = [(float(X[i, k]), float(loads[k])) for k in range(loads.size) if on[i, k]]
    return best_pos, best_m


def service_truck_loads(inp: Inputs) -> list[tuple[float, float]]:
    pos, _ = truck_position_max_midspan(inp)
    s = truck_fraction(inp) * (1.0 + impact_factor(inp)) * inp.live_load_factor
    return [(x, s * P) for x, P in pos]


def midspan_moment_statics(L: float, w: float, points) -> float:
    il = lambda x: x / 2 if x <= L / 2 else (L - x) / 2
    return w * L * L / 8.0 + sum(P * il(x) for x, P in points)
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_loads.py -v`
Expected: 6 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/loads.py tests/test_loads.py src/crackfem/params.py docs/sources.md docs/DECISIONS.md
git commit -m "feat: staged dead load, HS20-44 x 1.30 truck placement, impact and distribution" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 13: `nonlinear.py` (direct iteration, Ch11)

**Files:**
- Create: `src/crackfem/nonlinear.py`
- Test: `tests/test_nonlinear.py`

- [ ] **Step 1: Write the failing test**

```python
import numpy as np
import pytest

from crackfem.model import BeamModel, simply_supported, with_loads
from crackfem.nonlinear import solve_secant
from crackfem.solve import midpoint_moments, moment_at, solve

L, EI0, W = 24.5, 7.5e9, 28e3
MREF = W * L**2 / 8


def softening(M):
    return EI0 / (1.0 + np.abs(M) / MREF)


def test_determinate_beam_converges_in_one_pass():
    m = with_loads(simply_supported(L, 20), q=-W)
    res = solve_secant(m, softening, EI0)
    assert res.iterations == 1
    assert res.history[0] < 1e-12


def test_constant_stiffness_reproduces_the_linear_solution():
    m = with_loads(simply_supported(L, 20), q=-W)
    res = solve_secant(m, lambda M: np.full(M.size, EI0), EI0)
    np.testing.assert_allclose(res.solution.D, solve(m, EI0).D, rtol=1e-12)


def test_indeterminate_beam_redistributes_moment_and_reaches_a_fixed_point():
    n = 20
    m = with_loads(BeamModel(x=np.linspace(0, L, n + 1), fixed_dofs=(0, 1, 2 * n)), q=-W)
    res = solve_secant(m, softening, EI0, tol=1e-10)
    assert res.iterations >= 2
    np.testing.assert_allclose(softening(res.M_mid), res.EI, rtol=1e-8)
    assert res.solution.R.sum() == pytest.approx(W * L, rel=1e-10)
    linear_fixed_end = moment_at(m, solve(m, EI0), 0, 0.0)
    assert abs(moment_at(m, res.solution, 0, 0.0) / linear_fixed_end - 1.0) > 1e-3
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_nonlinear.py -v`
Expected: FAIL — `No module named 'crackfem.nonlinear'`

- [ ] **Step 3: Write `src/crackfem/nonlinear.py`**

```python
"""Direct iteration with secant stiffness (course Ch11): solve, update each element's EI
from its own moment, solve again, until the moments stop changing."""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .model import BeamModel
from .solve import Solution, midpoint_moments, solve


@dataclass(frozen=True)
class NonlinearResult:
    solution: Solution
    EI: np.ndarray          # stiffness used in the final solve
    M_mid: np.ndarray       # element midpoint moments of the final solve
    iterations: int         # number of stiffness updates
    history: tuple          # relative moment change after each update


def solve_secant(model: BeamModel, stiffness_fn: Callable[[np.ndarray], np.ndarray], EI_start,
                 tol: float = 1e-6, max_iter: int = 100) -> NonlinearResult:
    EI = np.broadcast_to(np.asarray(EI_start, dtype=float), (model.n_el,)).copy()
    sol = solve(model, EI)
    M = midpoint_moments(model, sol)
    history = []
    for it in range(1, max_iter + 1):
        EI = np.asarray(stiffness_fn(M), dtype=float)
        sol = solve(model, EI)
        M_new = midpoint_moments(model, sol)
        change = float(np.max(np.abs(M_new - M)) / max(float(np.max(np.abs(M_new))), 1e-12))
        history.append(change)
        M = M_new
        if change <= tol:
            return NonlinearResult(sol, EI, M, it, tuple(history))
    raise RuntimeError(f"secant iteration did not converge in {max_iter} updates: {history[-3:]}")
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_nonlinear.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/nonlinear.py tests/test_nonlinear.py
git commit -m "feat: direct secant iteration; one pass on determinate beams" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 14: `rootfind.py` (Illinois bracketed secant)

**Files:**
- Create: `src/crackfem/rootfind.py`
- Test: `tests/test_rootfind.py`

- [ ] **Step 1: Write the failing test**

```python
import math

import pytest
from scipy.optimize import brentq

from crackfem.rootfind import illinois

CASES = [
    (lambda x: x**3 - 2 * x - 5, 2.0, 3.0),
    (lambda x: math.cos(x) - x, 0.0, 1.0),
    (lambda x: 2 * max(0.0, 0.6 - x) - 0.1, 0.0, 1.0),   # flat beyond 0.6, like w_pred(eta)
]


@pytest.mark.parametrize("f,a,b", CASES)
def test_matches_brent(f, a, b):
    assert illinois(f, a, b, xtol=1e-12).x == pytest.approx(brentq(f, a, b, xtol=1e-14), abs=1e-10)


def test_refuses_an_unbracketed_interval():
    with pytest.raises(ValueError):
        illinois(lambda x: x * x + 1, -1.0, 1.0)


def test_records_history_and_counts_iterations():
    r = illinois(lambda x: x - 0.3, 0.0, 1.0)
    assert r.iterations >= 1
    assert len(r.history) == r.iterations + 2
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_rootfind.py -v`
Expected: FAIL — `No module named 'crackfem.rootfind'`

- [ ] **Step 3: Write `src/crackfem/rootfind.py`**

```python
"""Illinois method: a secant step that always keeps the root bracketed (the 'safeguarded
secant' of the spec). When the same end survives twice, its f value is halved so that
both ends keep moving."""
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class RootResult:
    x: float
    iterations: int
    history: tuple      # (x, f(x)) for both starting ends and every iterate


def illinois(f: Callable[[float], float], a: float, b: float,
             xtol: float = 1e-10, max_iter: int = 300) -> RootResult:
    fa, fb = f(a), f(b)
    history = [(a, fa), (b, fb)]
    if fa == 0.0:
        return RootResult(a, 0, tuple(history))
    if fb == 0.0:
        return RootResult(b, 0, tuple(history))
    if fa * fb > 0.0:
        raise ValueError("root is not bracketed: f(a) and f(b) have the same sign")
    side = 0
    for it in range(1, max_iter + 1):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        history.append((c, fc))
        if fc == 0.0:
            return RootResult(c, it, tuple(history))
        if fc * fb > 0.0:          # root in [a, c]: c replaces b, a survives
            b, fb = c, fc
            if side == -1:
                fa *= 0.5
            side = -1
        else:                      # root in [c, b]: c replaces a, b survives
            a, fa = c, fc
            if side == +1:
                fb *= 0.5
            side = +1
        if abs(b - a) <= xtol:
            return RootResult(0.5 * (a + b), it, tuple(history))
    raise RuntimeError(f"Illinois did not converge in {max_iter} iterations")
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_rootfind.py -v`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/rootfind.py tests/test_rootfind.py
git commit -m "feat: Illinois bracketed secant root finder" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 15: `inverse.py`

**Files:**
- Create: `src/crackfem/inverse.py`
- Test: `tests/test_inverse.py`

- [ ] **Step 1: Write the failing test** (synthetic forward map with a known inverse: w = 1e-3·max(0, 0.6 − η) m)

```python
import numpy as np
import pytest

from crackfem.inverse import estimate_eta, eta_decompression, is_monotone_nonincreasing


def w_pred(eta):
    return 1e-3 * max(0.0, 0.6 - eta)


@pytest.mark.parametrize("eta_true", [0.15, 0.4, 0.55])
def test_round_trip(eta_true):
    est = estimate_eta(w_pred, w_pred(eta_true))
    assert est.status == "ok"
    assert est.eta == pytest.approx(eta_true, abs=1e-6)


def test_no_crack_returns_a_bound_not_a_value():
    est = estimate_eta(w_pred, 0.0)
    assert est.status == "no_crack_bound"
    assert est.eta is None
    assert est.eta_lower_bound == pytest.approx(0.6, abs=1e-5)


def test_crack_wider_than_the_model_can_explain():
    est = estimate_eta(w_pred, 1.0e-3)
    assert est.status == "exceeds_model"
    assert est.eta is None


def test_intact_girder_that_already_cracks_is_reported_as_inconsistent():
    with pytest.raises(ValueError):
        estimate_eta(lambda e: 1e-3 * (1.2 - e), 0.1e-3)


def test_decompression_eta_and_monotonicity():
    assert eta_decompression(w_pred) == pytest.approx(0.6, abs=1e-5)
    assert is_monotone_nonincreasing(w_pred, np.linspace(0.02, 1.0, 50))
    assert not is_monotone_nonincreasing(lambda e: e, np.linspace(0.02, 1.0, 50))
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_inverse.py -v`
Expected: FAIL — `No module named 'crackfem.inverse'`

- [ ] **Step 3: Write `src/crackfem/inverse.py`**

```python
"""Inverse problem: remaining prestress fraction eta from a measured crack width.

w_pred(eta) is the forward model (metres); it must not increase with eta.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .rootfind import illinois


@dataclass(frozen=True)
class EtaEstimate:
    status: str                    # "ok" | "no_crack_bound" | "exceeds_model"
    eta: float | None
    eta_lower_bound: float | None  # set for "no_crack_bound": eta >= this value
    iterations: int


def eta_decompression(w_pred: Callable[[float], float], eta_lo: float = 0.02, eta_hi: float = 1.0,
                      w_tiny: float = 1e-9, xtol: float = 1e-7) -> float:
    """Smallest eta at which the predicted crack width is (numerically) zero."""
    g = lambda e: w_pred(e) - w_tiny
    if g(eta_hi) > 0.0:
        raise ValueError("the intact girder (eta_hi) already predicts a crack")
    if g(eta_lo) <= 0.0:
        return eta_lo
    return illinois(g, eta_lo, eta_hi, xtol=xtol).x


def estimate_eta(w_pred: Callable[[float], float], w_meas: float, eta_lo: float = 0.02,
                 eta_hi: float = 1.0, xtol: float = 1e-7) -> EtaEstimate:
    if w_meas <= 0.0:
        return EtaEstimate("no_crack_bound", None, eta_decompression(w_pred, eta_lo, eta_hi), 0)
    f = lambda e: w_pred(e) - w_meas
    if f(eta_hi) > 0.0:
        raise ValueError("the intact girder predicts a wider crack than measured: model inconsistent")
    if f(eta_lo) < 0.0:
        return EtaEstimate("exceeds_model", None, None, 0)
    r = illinois(f, eta_lo, eta_hi, xtol=xtol)
    return EtaEstimate("ok", r.x, None, r.iterations)


def is_monotone_nonincreasing(w_pred: Callable[[float], float], grid) -> bool:
    w = np.array([w_pred(float(e)) for e in grid])
    return bool(np.all(np.diff(w) <= 1e-15))
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_inverse.py -v`
Expected: 7 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/inverse.py tests/test_inverse.py
git commit -m "feat: eta inverse with decompression bound and model-range checks" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 15b: Staged construction — slab locked-in strain + `staging.py`

**Files:**
- Modify: `src/crackfem/section.py`
- Create: `src/crackfem/staging.py`
- Test: `tests/test_staging.py`

Decision 2026-10-12 (student): model the precast girder in two stages. **Stage A:** the bare girder carries the prestress and the weight of the girder and the wet slab. **Stage B:** after the slab hardens, the wearing surface, barriers and trucks act on the composite section. The slab is cast stress-free on the deformed girder, so each slab fibre has a *stress-free strain* equal to the stage-A strain plane at its level. A composite solve at `M = M_A` then reproduces stage A with an unstressed slab, and moments above `M_A` are stage-B increments on the composite section. Concrete stress everywhere becomes `E · (ε(y) − ε_free)`.

- [ ] **Step 1: Write the failing test** — `tests/test_staging.py` (independent evidence = the classic two-stage hand calculation with transformed sections built from `gross_properties`, not from fibres)

```python
import numpy as np
import pytest

from crackfem.geometry import discretize, gross_properties, slab
from crackfem.section import (PCSection, concrete_stress, resultants, soffit_stress, solve_state,
                              uncracked_EI)
from crackfem.staging import design_prestress_staged, girder_only, staged_section

E_G, E_S, E_P = 30e9, 27e9, 195e9
H_G, Y_P, A_P = 0.8, 0.1, 1.5e-3
GIRDER = slab(0.0, H_G, 0.4, E=E_G, fct=3e6)        # a rectangle is a trapezoid with b0 == b1
DECK = slab(H_G, 0.15, 1.2, E=E_S, fct=2.5e6)
FIB = discretize([GIRDER, DECK], 2e-4)
IS_GIRDER = FIB.y < H_G
H = H_G + 0.15
E0, M_A, DM = 4e-3, 200e3, 300e3


def composite(A_p=A_P):
    return PCSection(FIB, Y_P, A_p, E_P, H)


def transformed(parts):
    """Transformed section in girder concrete (n = E / E_G), strand included."""
    A, y, I = gross_properties(parts, E_ref=E_G)
    n = E_P / E_G
    At = A + n * A_P
    yt = (A * y + n * A_P * Y_P) / At
    It = I + A * (y - yt) ** 2 + n * A_P * (Y_P - yt) ** 2
    return At, yt, It


def test_composite_solve_at_the_stage_a_moment_reproduces_stage_a_with_an_unstressed_slab():
    sec = staged_section(composite(), IS_GIRDER, H_G, E0, M_A)
    a = solve_state(girder_only(composite(), IS_GIRDER, H_G), E0, M_A, cracked=False)
    b = solve_state(sec, E0, M_A, cracked=False)
    assert (b.eps0, b.kappa) == pytest.approx((a.eps0, a.kappa), rel=1e-9)
    assert np.abs(concrete_stress(sec, b)[~IS_GIRDER]).max() < 1.0      # Pa


def test_stage_b_increment_matches_the_two_stage_hand_calculation():
    At_g, yt_g, It_g = transformed([GIRDER])
    At_c, yt_c, It_c = transformed([GIRDER, DECK])
    F0 = E_P * A_P * E0
    e_g = yt_g - Y_P
    soffit_a = -F0 / At_g - F0 * e_g * yt_g / It_g + M_A * yt_g / It_g     # stage A, bare girder
    soffit_b = DM * yt_c / It_c                                          # stage B, composite
    sec = staged_section(composite(), IS_GIRDER, H_G, E0, M_A)
    a = solve_state(sec, E0, M_A, cracked=False)
    b = solve_state(sec, E0, M_A + DM, cracked=False)
    assert soffit_stress(sec, b) == pytest.approx(soffit_a + soffit_b, rel=1e-5)
    assert b.kappa - a.kappa == pytest.approx(DM / (E_G * It_c), rel=1e-5)
    assert uncracked_EI(sec) == pytest.approx(E_G * It_c, rel=1e-5)


@pytest.mark.parametrize("M", [1.2e6, 1.8e6])
def test_cracked_staged_state_is_in_equilibrium_with_no_concrete_tension(M):
    sec = staged_section(composite(), IS_GIRDER, H_G, E0, M_A)
    st = solve_state(sec, E0, M, cracked=True)
    N, M_int = resultants(sec, E0, st)
    assert abs(N) < 1e-6 * E_P * A_P * E0
    assert M_int == pytest.approx(M, rel=1e-6)
    assert concrete_stress(sec, st).max() <= 0.0


def test_staged_design_is_the_smallest_strand_count_with_no_soffit_tension():
    d = design_prestress_staged(composite, IS_GIRDER, H_G, 98.7e-6, 1.0e9, M_A, M_A + DM)
    assert d.soffit_stress_service <= 0.0
    fewer = composite((d.n_strands - 1) * 98.7e-6)
    from crackfem.prestress import eps_p0_for_stress
    e0 = eps_p0_for_stress(girder_only(fewer, IS_GIRDER, H_G), 1.0e9)
    sec = staged_section(fewer, IS_GIRDER, H_G, e0, M_A)
    assert soffit_stress(sec, solve_state(sec, e0, M_A + DM, cracked=False)) > 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_staging.py -v`
Expected: FAIL — `No module named 'crackfem.staging'`

- [ ] **Step 3: Add the stress-free strain to `src/crackfem/section.py`.** Make exactly these changes and keep everything else (docstrings, the `State` class, `strand_stress`, `delta_sigma_p`, `resultants`, `uncracked_EI`, `cracking_moment`) as it is:

(a) Add one sentence at the end of the module docstring: `Staged sections: each fibre may carry a stress-free strain eps_free (the slab's strain locked in at casting); concrete stress is E * (eps(y) - eps_free).`

(b) Give `PCSection` a last field and a helper method (keep its existing `soffit_E` property; if it reads the soffit E by `np.argmin(self.fib.y)`, keep that):

```python
    eps_free: np.ndarray | None = None      # per-fibre stress-free strain; None = all zero

    def free_strain(self) -> np.ndarray:
        return np.zeros(self.fib.y.size) if self.eps_free is None else self.eps_free
```

(c) Replace `_solve_active` with:

```python
def _solve_active(sec: PCSection, eps_p0: float, M: float, active: np.ndarray) -> tuple[float, float]:
    f = sec.fib
    ea = f.E[active] * f.A[active]
    y = f.y[active]
    ef = sec.free_strain()[active]
    epap = sec.E_p * sec.A_p
    EA = ea.sum() + epap
    ES = (ea * y).sum() + epap * sec.y_p
    EI = (ea * y * y).sum() + epap * sec.y_p**2
    J = np.array([[EA, -ES], [-ES, EI]])
    rhs = np.array([-epap * eps_p0 + (ea * ef).sum(),
                    M + epap * sec.y_p * eps_p0 - (ea * ef * y).sum()])
    eps0, kappa = np.linalg.solve(J, rhs)
    return float(eps0), float(kappa)
```

(d) In `solve_state`: change the trivial-state guard to `if eps_p0 == 0.0 and M == 0.0 and sec.eps_free is None:`, and change the active-set test to use the mechanical strain:

```python
        new_active = (eps0 - kappa * sec.fib.y - sec.free_strain()) < 0.0
```

(e) Replace `concrete_stress` with:

```python
def concrete_stress(sec: PCSection, st: State) -> np.ndarray:
    sig = sec.fib.E * (st.strain(sec.fib.y) - sec.free_strain())
    return np.minimum(sig, 0.0) if st.cracked else sig
```

(f) In `soffit_stress`, keep the cracked-state guard and return `sec.soffit_E * (st.eps0 - float(sec.free_strain()[int(np.argmin(sec.fib.y))]))`.

- [ ] **Step 4: Run the existing section tests (they must still pass unchanged)**

Run: `python -m pytest tests/test_section_uncracked.py tests/test_section_cracked.py tests/test_prestress.py -q`
Expected: all passed (21)

- [ ] **Step 5: Write `src/crackfem/staging.py`**

```python
"""Staged construction of a precast prestressed girder with a cast-in-place slab (unshored).

Stage A: the bare girder carries the prestress and the weight of the girder and the wet slab.
Stage B: the hardened slab acts compositely; later loads act on the composite section.
The slab is cast stress-free on the deformed girder, so its stress-free strain is the
stage-A strain plane at its level. A composite solve at M = M_A therefore reproduces
stage A with an unstressed slab, and moments above M_A are stage-B increments.
"""
from dataclasses import replace
from typing import Callable

import numpy as np

from .geometry import Fibers
from .prestress import PrestressDesign, eps_p0_for_stress
from .section import PCSection, soffit_stress, solve_state


def girder_only(sec: PCSection, is_girder: np.ndarray, h_girder: float) -> PCSection:
    """The bare girder: the girder fibres of a composite section, same strand."""
    f = sec.fib
    fib = Fibers(f.y[is_girder], f.A[is_girder], f.E[is_girder], f.fct[is_girder])
    return PCSection(fib, sec.y_p, sec.A_p, sec.E_p, h_girder)


def staged_section(sec: PCSection, is_girder: np.ndarray, h_girder: float,
                   eps_p0_construction: float, M_A: float) -> PCSection:
    """Composite section whose slab fibres carry the stage-A strain plane as stress-free strain."""
    plane = solve_state(girder_only(sec, is_girder, h_girder), eps_p0_construction, M_A, cracked=False)
    return replace(sec, eps_free=np.where(is_girder, 0.0, plane.strain(sec.fib.y)))


def design_prestress_staged(make_section: Callable[[float], PCSection], is_girder: np.ndarray,
                            h_girder: float, strand_area: float, fpe: float, M_A: float,
                            M_total: float, n_max: int = 150) -> PrestressDesign:
    """Smallest strand count with no soffit tension under the total service moment, with stage A
    on the bare girder. fpe = effective strand stress on the bare girder at zero external moment."""
    for n in range(1, n_max + 1):
        comp = make_section(n * strand_area)
        e0 = eps_p0_for_stress(girder_only(comp, is_girder, h_girder), fpe)
        sec = staged_section(comp, is_girder, h_girder, e0, M_A)
        s = soffit_stress(sec, solve_state(sec, e0, M_total, cracked=False))
        if s <= 0.0:
            return PrestressDesign(n, n * strand_area, fpe, e0, fpe * n * strand_area, M_total, s)
    raise ValueError(f"no strand count up to {n_max} removes soffit tension")
```

- [ ] **Step 6: Run to verify it passes, then the full suite**

Run: `python -m pytest tests/test_staging.py -v` → Expected: 5 passed
Run: `python -m pytest -q` → Expected: all passed

If `test_stage_b_increment_matches_the_two_stage_hand_calculation` fails, the sign of the `eps_free` terms in `_solve_active` is the first suspect — re-derive N = Σ E A (ε0 − κy − ε_free) + Ep Ap (ε_p0 + ε0 − κ y_p) = 0 and M = −Σ σ A y; do not edit the test.

- [ ] **Step 7: Commit**

```bash
git add src/crackfem/section.py src/crackfem/staging.py tests/test_staging.py
git commit -m "feat: staged construction (slab locked-in strain, staged prestress design)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---


## Task 16: `girder1d.py` (Stage-1 wiring, staged construction)

**Files:**
- Create: `src/crackfem/girder1d.py`
- Test: `tests/test_girder1d.py`

Depends on Tasks 9 (`ec2_secant_EI`/`branson_secant_EI` with `M_ref`, `ec2_kappa`), 11 (`node_moments`), 12 (`stage_a_line`, `stage_b_line`, live-load factor) and 15b (`staging.py`). Stage A moments come from statics (the bare girder is determinate); the FE model carries the stage-B loads (wearing surface + barriers, then + truck) on the composite girder. Each element and each node gets its own staged section because the slab's locked-in strain depends on the local stage-A moment. η scales today's strand strain; the slab strain locked in at construction (η = 1) does not change.

- [ ] **Step 1: Write the failing test**

```python
import numpy as np
import pytest

from crackfem.girder1d import Girder1D
from crackfem.inverse import estimate_eta, eta_decompression, is_monotone_nonincreasing
from crackfem.solve import solve


@pytest.fixture(scope="module")
def g():
    return Girder1D()


def test_total_dead_moment_is_stage_a_plus_stage_b(g):
    L = g.inp.span
    assert g.M_dead_mid == pytest.approx((g.q_a + g.q_b) * L**2 / 8, rel=1e-9)


def test_intact_girder_shows_no_crack_at_rest(g):
    assert g.w_pred(1.0) == 0.0


def test_predicted_width_does_not_increase_with_eta(g):
    assert is_monotone_nonincreasing(g.w_pred, np.linspace(0.02, 1.0, 40))


def test_round_trip_through_the_real_section(g):
    eta_true = round(eta_decompression(g.w_pred) - 0.15, 3)
    est = estimate_eta(g.w_pred, g.w_pred(eta_true))
    assert est.status == "ok"
    assert est.eta == pytest.approx(eta_true, abs=1e-5)


def test_intact_service_response_is_the_linear_uncracked_solution(g):
    ind = g.indicator(1.0)
    lin = solve(g.model_service, g.EI_I)
    assert ind["deflection_service_mm"] == pytest.approx(-lin.D[2 * g.mid_node] * 1e3, rel=1e-9)
    assert all(r == pytest.approx(1.0) for r in ind["EI_ratio"])
    assert ind["visible_crack_zone_m"] == 0.0


def test_damaged_girder_is_softer_but_still_one_pass(g):
    eta = round(eta_decompression(g.w_pred) - 0.15, 3)
    ind, ref = g.indicator(eta), g.indicator(1.0)
    assert ind["deflection_service_mm"] > ref["deflection_service_mm"]
    assert any(ind["precracked"])
    assert min(ind["EI_ratio"]) < 1.0
    assert ind["iterations_service"] == 1
    assert 0.0 < ind["visible_crack_zone_m"] < g.inp.span


def test_camber_is_upward_and_stage_a_deflection_downward(g):
    ind = g.indicator(1.0)
    assert ind["camber_mm"] > 0.0
    assert ind["deflection_stage_a_mm"] > 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_girder1d.py -v`
Expected: FAIL — `No module named 'crackfem.girder1d'`

- [ ] **Step 3: Write `src/crackfem/girder1d.py`**

```python
"""Stage 1: the whole chain on one simply supported girder (Gk), built in two stages.

Stage A: the bare girder carries the prestress and the weight of the girder and the wet slab
(moments from statics). Stage B: wearing surface, barriers and the truck act on the composite
girder through the FE model. Each element and node has its own staged section, because the
slab's locked-in strain depends on the local stage-A moment. eta scales today's strand strain;
the slab strain locked in at construction (eta = 1) does not change.

Crack width at rest (total dead load) -> eta through the inverse; stiffness per element from
its own total moment, measured from the end of stage A; deflection under stage-B dead load +
truck. Elements whose total service moment has exceeded Mcr(eta) are 'pre-cracked' (traffic
cracked them): their tensile strength is taken as zero, so their effective Mcr is the
decompression moment.
"""
import numpy as np

from .crackwidth import EC2Crack, effective_tension_height
from .geometry import area_below, discretize, girder_outline, gross_properties, slab, top_of
from .loads import midspan_moment_statics, service_truck_loads, stage_a_line, stage_b_line
from .materials import ec2_ecm, ec2_fctm
from .model import simply_supported, with_loads
from .nonlinear import solve_secant
from .params import Inputs
from .prestress import effective_prestress
from .section import PCSection, cracking_moment, delta_sigma_p, solve_state, uncracked_EI
from .solve import midpoint_moments, node_moments, solve
from .staging import design_prestress_staged, girder_only, staged_section
from .stiffness import branson_secant_EI, ec2_kappa, ec2_secant_EI


class Girder1D:
    def __init__(self, inp: Inputs = Inputs(), n_el: int = 40, dy_max: float = 1e-3):
        if n_el % 2:
            raise ValueError("n_el must be even so that midspan is a node")
        self.inp, self.n_el = inp, n_el
        L = inp.span
        self.Ec_girder = ec2_ecm(inp.fck_girder)
        self.fct_girder = ec2_fctm(inp.fck_girder)
        self.girder_parts = girder_outline(inp.girder_outline_in, self.Ec_girder, self.fct_girder)
        self.h_girder = top_of(self.girder_parts)
        self.parts = self.girder_parts + [slab(self.h_girder, inp.slab_thickness, inp.girder_spacing,
                                               ec2_ecm(inp.fck_slab), ec2_fctm(inp.fck_slab))]
        self.h = top_of(self.parts)
        self.fibers = discretize(self.parts, dy_max)
        self.is_girder = self.fibers.y < self.h_girder
        self.q_a = stage_a_line(inp, gross_properties(self.girder_parts)[0])
        self.q_b = stage_b_line(inp)
        self.truck = service_truck_loads(inp)
        base = simply_supported(L, n_el)
        self.mid_node = n_el // 2

        def stage_a(x):                     # statics: the bare girder is simply supported
            return self.q_a * x * (L - x) / 2.0

        self.M_A_el, self.M_A_nd = stage_a(base.midpoints), stage_a(base.x)
        self.M_service = midspan_moment_statics(L, self.q_a + self.q_b, self.truck)
        fpe = effective_prestress(inp.fpj_ratio, inp.fpu, inp.lump_sum_loss)
        self.design = design_prestress_staged(self.make_section, self.is_girder, self.h_girder,
                                              inp.strand_area, fpe, float(self.M_A_nd[self.mid_node]),
                                              self.M_service)
        self.composite = self.make_section(self.design.A_p)
        e_c = self.design.eps_p0            # construction state: eta = 1
        self.sec_el = [staged_section(self.composite, self.is_girder, self.h_girder, e_c, float(m))
                       for m in self.M_A_el]
        self.sec_nd = [staged_section(self.composite, self.is_girder, self.h_girder, e_c, float(m))
                       for m in self.M_A_nd]
        self.EI_I = uncracked_EI(self.composite)
        self.EI_girder = uncracked_EI(girder_only(self.composite, self.is_girder, self.h_girder))
        self.crack = EC2Crack(inp.k1, inp.k2, inp.k3, inp.k4, inp.kt, inp.bond_xi, inp.Ep,
                              self.Ec_girder, self.fct_girder, inp.strand_cover,
                              inp.strand_diameter, self.design.A_p)
        self.model_dead = with_loads(base, q=-self.q_b)
        self.model_service = with_loads(base, q=-self.q_b,
                                        point_loads=[(x, -P) for x, P in self.truck])
        # Determinate girder: stage-B moments do not depend on stiffness (test_nonlinear proves the
        # secant loop stops after one pass), so the total moments are computed once.
        self.M_dead_nd = self.M_A_nd + node_moments(self.model_dead, solve(self.model_dead, self.EI_I))
        self.M_dead_mid = float(self.M_dead_nd[self.mid_node])
        self.M_service_el = self.M_A_el + midpoint_moments(self.model_service,
                                                           solve(self.model_service, self.EI_I))

    def make_section(self, A_p: float) -> PCSection:
        return PCSection(self.fibers, self.inp.strand_centroid, A_p, self.inp.Ep, self.h)

    def eps_p0(self, eta: float) -> float:
        return eta * self.design.eps_p0

    def mcr(self, eta: float, sec: PCSection, fct: float | None = None) -> float:
        """Total moment at which this staged section's soffit reaches fct (0 = decompression)."""
        return cracking_moment(sec, self.eps_p0(eta), self.fct_girder if fct is None else fct)

    def crack_at(self, eta: float, M: float, sec: PCSection) -> dict:
        st = solve_state(sec, self.eps_p0(eta), M, cracked=True)
        ds = delta_sigma_p(sec, st)
        if ds <= 0.0:
            return {"delta_sigma_p": ds, "Ac_eff": None, "width": 0.0, "neutral_axis": None}
        hce = effective_tension_height(self.h, sec.y_p, st.neutral_axis)
        Ac = area_below(self.girder_parts, hce)
        return {"delta_sigma_p": ds, "Ac_eff": Ac, "width": self.crack.width(ds, Ac),
                "neutral_axis": st.neutral_axis}

    def w_pred(self, eta: float) -> float:
        """Predicted crack width (m) at midspan under the total dead load (the inspection state)."""
        return self.crack_at(eta, self.M_dead_mid, self.sec_nd[self.mid_node])["width"]

    def precracked_zone(self, eta: float) -> np.ndarray:
        return np.array([m > self.mcr(eta, s) for m, s in zip(self.M_service_el, self.sec_el)])

    def stiffness_fn(self, eta: float, precracked: np.ndarray, method: str = "ec2"):
        """EI per element from its stage-B moment (what the FE model returns); the element's
        total moment is M_A + M_B and its reference state is the end of stage A."""
        e0 = self.eps_p0(eta)
        mcr = [self.mcr(eta, s, 0.0 if pre else None) for s, pre in zip(self.sec_el, precracked)]

        def fn(M_b: np.ndarray) -> np.ndarray:
            out = np.empty(M_b.size)
            for i, mb in enumerate(M_b):
                M_ref = float(self.M_A_el[i])
                M = M_ref + float(mb)
                if method == "ec2":
                    out[i] = ec2_secant_EI(self.sec_el[i], e0, M, mcr[i], self.inp.beta, M_ref)
                else:
                    out[i] = branson_secant_EI(self.sec_el[i], e0, M, mcr[i], M_ref)
            return out

        return fn

    def analyse(self, eta: float, model, method: str = "ec2"):
        return solve_secant(model, self.stiffness_fn(eta, self.precracked_zone(eta), method), self.EI_I)

    def kappa_curves(self, eta: float, M_grid) -> dict:
        """Midspan moment-curvature for the report figure: total moment, curvature measured from
        the uncracked state at the end of stage A. The fully cracked curve is drawn only above
        the decompression moment (below it, prestress alone would crack the top)."""
        sec = self.sec_nd[self.mid_node]
        e0 = self.eps_p0(eta)
        k_ref = solve_state(sec, e0, float(self.M_A_nd[self.mid_node]), cracked=False).kappa
        mcr_sound, mcr_pre = self.mcr(eta, sec), self.mcr(eta, sec, 0.0)
        out = {"uncracked": [], "cracked": [], "ec2_sound": [], "ec2_precracked": []}
        for M in M_grid:
            M = float(M)
            out["uncracked"].append(solve_state(sec, e0, M, cracked=False).kappa - k_ref)
            out["cracked"].append(solve_state(sec, e0, M, cracked=True).kappa - k_ref
                                  if M >= mcr_pre else float("nan"))
            out["ec2_sound"].append(ec2_kappa(sec, e0, M, mcr_sound, self.inp.beta) - k_ref)
            out["ec2_precracked"].append(ec2_kappa(sec, e0, M, mcr_pre, self.inp.beta) - k_ref)
        return out

    def visible_crack_zone(self, eta: float) -> float:
        """Length of span (m) whose nodes show an open crack at rest. On synthetic data this is a
        model output, not a held-out check: the masks come from the same model."""
        open_ = [self.crack_at(eta, float(M), s)["width"] > 0.0 for M, s in zip(self.M_dead_nd, self.sec_nd)]
        xs = self.model_dead.x[open_]
        return float(xs.max() - xs.min()) if xs.size else 0.0

    def indicator(self, eta: float, method: str = "ec2") -> dict:
        dead = self.analyse(eta, self.model_dead, method)
        serv = self.analyse(eta, self.model_service, method)
        v_dead = -dead.solution.D[2 * self.mid_node]
        v_serv = -serv.solution.D[2 * self.mid_node]
        L = self.inp.span
        sec_mid = self.sec_nd[self.mid_node]
        k_transfer = solve_state(girder_only(self.composite, self.is_girder, self.h_girder),
                                 self.design.eps_p0, 0.0, cracked=False).kappa
        return {
            "eta": eta,
            "method": method,
            "w_pred_mm": self.w_pred(eta) * 1e3,
            "Mcr_mid_kNm": self.mcr(eta, sec_mid) / 1e3,
            "Mdec_mid_kNm": self.mcr(eta, sec_mid, 0.0) / 1e3,
            "M_dead_mid_kNm": self.M_dead_mid / 1e3,
            "camber_mm": -k_transfer * L * L / 8.0 * 1e3,            # bare girder at transfer, eta = 1
            "deflection_stage_a_mm": 5.0 * self.q_a * L**4 / (384.0 * self.EI_girder) * 1e3,
            "deflection_dead_mm": v_dead * 1e3,                      # stage-B dead load, composite
            "deflection_service_mm": v_serv * 1e3,                   # stage-B dead load + truck
            "deflection_truck_mm": (v_serv - v_dead) * 1e3,
            "deflection_limit_mm": L / self.inp.deflection_limit_ratio * 1e3,
            "visible_crack_zone_m": self.visible_crack_zone(eta),
            "EI_ratio": (serv.EI / self.EI_I).tolist(),
            "x_mid": self.model_service.midpoints.tolist(),
            "x_nodes": self.model_service.x.tolist(),
            "v_service_mm": (-serv.solution.D[0::2] * 1e3).tolist(),
            "precracked": self.precracked_zone(eta).tolist(),
            "iterations_service": serv.iterations,
        }
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_girder1d.py -v`
Expected: 7 passed. If `test_intact_girder_shows_no_crack_at_rest` fails, check `design.soffit_stress_service <= 0` first, then the sign of `delta_sigma_p`.

- [ ] **Step 5: Run the whole suite**

Run: `python -m pytest -q`
Expected: all passed

- [ ] **Step 6: Progress log entry B2** — append to `docs/PROGRESS_LOG.md` (fill every line from the actual test run; design numbers such as strand count are reported in Entry 3 after `run_all.py`):

```markdown
## Entry 2 — B2 Stage-1 core (<date>)

**Planned:** section, crack width, stiffness, beam FE, nonlinear loop, inverse — each with independent tests (spec §6).
**Found:** <test findings worth reporting, e.g. "beam nodal deflections exact for N = 1…40 (proposal P2 confirmed)", "secant loop converges in 1 pass on the determinate girder (P5 confirmed)", "indeterminate toy beam needs N passes", "a section whose strand is below the kern cracks at the top under prestress alone (test corrected)", "staged section matches the two-stage hand calculation">.
**Changed & why:** staged construction adopted (student decision 2026-10-12, reviewer estimate ≈ 0.15 shift in η if ignored); HS20-44 × 1.30 applied (student decision 2026-10-12); <any other deviation, with reason>.
**Numbers:** `python -m pytest -q` → <count> passed, <count> skipped.
**Reading guide:** read in this order: elements.py → model.py → solve.py → nonlinear.py (the course part), then section.py → staging.py → crackwidth.py → stiffness.py → girder1d.py (the physics). Explain-back: (1) why are the nodal deflections exact for any N? (2) why does the secant loop stop after one pass here? (3) why does the slab carry no stress from its own weight, and what does that change in the crack width at rest?
**Reviewed by student:** pending (after 2026-10-19).
```

(Replace every `<...>` with actual content; nothing in angle brackets may remain.)

- [ ] **Step 7: Commit**

```bash
git add src/crackfem/girder1d.py tests/test_girder1d.py docs/PROGRESS_LOG.md
git commit -m "feat: Stage-1 staged girder wiring (w_pred, stiffness, indicator)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---


## Task 17: `uq.py`

**Files:**
- Create: `src/crackfem/uq.py`
- Test: `tests/test_uq.py`

- [ ] **Step 1: Write the failing test** (independent check = for the linear map η = 0.6 − w/1e-3 the η percentiles are analytic)

```python
import numpy as np
import pytest

from crackfem.params import Inputs
from crackfem.uq import monte_carlo_eta, sample_widths, summarize, tornado


def w_pred(eta):
    return 1e-3 * max(0.0, 0.6 - eta)


def test_samples_are_positive_reproducible_and_sized():
    a = sample_widths(0.2e-3, 0.1e-3, 500, seed=1)
    b = sample_widths(0.2e-3, 0.1e-3, 500, seed=1)
    assert a.size == 500 and (a > 0).all()
    np.testing.assert_array_equal(a, b)


def test_zero_uncertainty_collapses_to_the_point_estimate():
    mc = monte_carlo_eta(w_pred, 0.2e-3, 0.0, n=20)
    assert np.allclose(mc.etas, 0.4, atol=1e-6)


def test_linear_map_gives_the_analytic_band():
    sigma_w = 0.1e-3 / 2
    mc = monte_carlo_eta(w_pred, 0.2e-3, 0.1e-3, n=20000, seed=7)
    s = mc.summary()
    assert s["p50"] == pytest.approx(0.4, abs=0.002)
    assert s["p95"] - s["p05"] == pytest.approx(2 * 1.6449 * sigma_w / 1e-3, rel=0.05)
    assert mc.status_counts == {"ok": 20000}


def test_summary_percentiles_are_ordered():
    s = summarize(np.random.default_rng(3).normal(size=1000))
    assert s["p05"] < s["p50"] < s["p95"]


def test_tornado_sorts_by_swing():
    metric = lambda inp: inp.girder_spacing * 10 + inp.slab_thickness
    rows = tornado(metric, Inputs(), {"slab_thickness": (0.18, 0.22), "girder_spacing": (2.0, 2.4)})
    assert [r["name"] for r in rows] == ["girder_spacing", "slab_thickness"]
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_uq.py -v`
Expected: FAIL — `No module named 'crackfem.uq'`

- [ ] **Step 3: Write `src/crackfem/uq.py`**

```python
"""Monte Carlo on the measured crack width, and one-at-a-time sensitivity (tornado).

ClaimBridge reports w +/- u with u = 2 GSD; u is read as a 95% bound, so sigma = u / 2.
Widths are truncated at w > 0 by rejection.
"""
from dataclasses import dataclass, replace
from typing import Callable

import numpy as np

from .inverse import estimate_eta
from .params import Inputs


def sample_widths(w_meas: float, u95: float, n: int, seed: int = 20261101) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(0)
    while out.size < n:
        draw = rng.normal(w_meas, u95 / 2.0, size=2 * n)
        out = np.concatenate([out, draw[draw > 0.0]])
    return out[:n]


def summarize(values) -> dict:
    v = np.asarray(values, dtype=float)
    p05, p50, p95 = np.percentile(v, [5, 50, 95])
    return {"n": int(v.size), "mean": float(v.mean()), "p05": float(p05), "p50": float(p50), "p95": float(p95)}


@dataclass(frozen=True)
class MonteCarlo:
    widths: np.ndarray
    etas: np.ndarray            # only samples with status "ok"
    status_counts: dict

    def summary(self) -> dict | None:
        return summarize(self.etas) if self.etas.size else None


def monte_carlo_eta(w_pred: Callable[[float], float], w_meas: float, u95: float, n: int = 2000,
                    seed: int = 20261101, **estimate_kw) -> MonteCarlo:
    widths = np.full(n, w_meas) if u95 == 0.0 else sample_widths(w_meas, u95, n, seed)
    etas, counts = [], {}
    for w in widths:
        est = estimate_eta(w_pred, float(w), **estimate_kw)
        counts[est.status] = counts.get(est.status, 0) + 1
        if est.status == "ok":
            etas.append(est.eta)
    return MonteCarlo(widths, np.array(etas), counts)


def tornado(metric: Callable[[Inputs], float], base: Inputs, perturbations: dict) -> list[dict]:
    """perturbations: {field name: (low value, high value)}; rows sorted by swing."""
    m0 = metric(base)
    rows = []
    for name, (lo, hi) in perturbations.items():
        rows.append({"name": name, "low": lo, "high": hi, "base": m0,
                     "metric_low": metric(replace(base, **{name: lo})),
                     "metric_high": metric(replace(base, **{name: hi}))})
    return sorted(rows, key=lambda r: -abs(r["metric_high"] - r["metric_low"]))
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_uq.py -v`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/uq.py tests/test_uq.py
git commit -m "feat: Monte Carlo on crack width and tornado sensitivity" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 18: `claimbridge.py` + live pipeline tests

**Files:**
- Create: `src/crackfem/claimbridge.py`, `tests/fixtures/spike_ok.txt`, `tests/fixtures/spike_refused.txt`
- Test: `tests/test_claimbridge.py`

- [ ] **Step 1: Write the fixtures.** `tests/fixtures/spike_ok.txt` is the real output captured on 2026-10-07 (paths shortened to file names):

```text
Mask:      mask.png
Surface:   s.json
GSD:       0.0500 mm/px  (min measurable width 0.15 mm)
--- crack 1: 394 pts (0 refused, 0 flagged)
    max width: 6.00 px = 0.30 ± 0.10 mm
    DOH rating (rule engine): 5 — ดีมาก (Excellent)
```

`tests/fixtures/spike_refused.txt` follows the pipeline's report format (`REFUSED: <gate>`):

```text
Mask:      mask.png
Surface:   s.json
GSD:       0.3800 mm/px  (min measurable width 1.14 mm)
--- crack 1: 390 pts (0 refused, 0 flagged)
    REFUSED: below-measurable-resolution
```

- [ ] **Step 2: Write the failing test**

```python
import json
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from crackfem.claimbridge import (default_config_path, governing_width, load_pipeline_config,
                                  parse_pipeline_output, run_pipeline, write_mask, write_surface)

FIX = Path(__file__).parent / "fixtures"


def test_mask_has_the_requested_crack_width(tmp_path):
    p = tmp_path / "m.png"
    write_mask(p, width_px=5, n_cracks=3)
    img = np.array(Image.open(p))
    row = img[img.shape[0] // 2]
    runs = np.diff(np.flatnonzero(np.diff(np.concatenate([[0], row > 127, [0]]).astype(int))))[::2]
    assert list(runs) == [5, 5, 5]


def test_surface_json_has_the_fields_the_pipeline_reads(tmp_path):
    p = tmp_path / "s.json"
    write_surface(p, 5e-5)
    meta = json.loads(p.read_text(encoding="utf-8"))
    assert meta["meters_per_pixel_x"] == meta["meters_per_pixel_y"] == 5e-5
    assert np.load(tmp_path / meta["mapping_outputs"]["xyz_npy"]).shape == (1, 1, 3)


def test_parse_real_output():
    out = parse_pipeline_output((FIX / "spike_ok.txt").read_text(encoding="utf-8"))
    assert out["gsd_mm"] == 0.05 and out["min_measurable_mm"] == 0.15
    c = out["cracks"][0]
    assert (c["gate"], c["max_width_px"], c["max_width_mm"], c["uncertainty_mm"]) == ("ok", 6.0, 0.30, 0.10)
    assert c["rating"] == "5" and c["rating_label"] == "ดีมาก (Excellent)"
    assert governing_width(out)["max_width_mm"] == 0.30


def test_parse_refused_output():
    out = parse_pipeline_output((FIX / "spike_refused.txt").read_text(encoding="utf-8"))
    assert out["cracks"][0]["gate"] == "below-measurable-resolution"
    assert governing_width(out) is None


@pytest.mark.skipif(not default_config_path().is_file(),
                    reason="needs a pipeline config (CRACK_PIPELINE_CONFIG or pipeline.local.json)")
def test_live_pipeline_measures_a_six_pixel_crack(tmp_path):
    write_mask(tmp_path / "m.png", width_px=6)
    write_surface(tmp_path / "s.json", 5e-5)
    out = parse_pipeline_output(run_pipeline(load_pipeline_config(), tmp_path / "m.png", tmp_path / "s.json"))
    assert governing_width(out)["max_width_mm"] == pytest.approx(0.30)
```

- [ ] **Step 3: Run to verify it fails**

Run: `python -m pytest tests/test_claimbridge.py -v`
Expected: FAIL — `No module named 'crackfem.claimbridge'`

- [ ] **Step 4: Write `src/crackfem/claimbridge.py`**

```python
"""Bridge to the author's thesis inspection pipeline (ClaimBridge, private).

Only files cross the boundary: a PNG mask + a surface JSON go in; the pipeline's text report
comes out and is parsed. How to run the pipeline is configured per machine in
pipeline.local.json (not committed; see pipeline.example.json), so this repository contains
no detail of the private code base.

The mask is a PNG with luminance > 127 = crack; the surface JSON fixes the GSD.
"""
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image


def write_mask(path, width_px: int, n_cracks: int = 3, length_px: int = 400,
               spacing_px: int = 150, margin_px: int = 50) -> None:
    """Vertical flexural cracks of constant width, white on black."""
    if width_px < 1:
        raise ValueError("width_px must be at least 1")
    H = length_px + 2 * margin_px
    W = 2 * margin_px + spacing_px * (n_cracks - 1) + width_px
    img = np.zeros((H, W), dtype=np.uint8)
    for i in range(n_cracks):
        x0 = margin_px + i * spacing_px
        img[margin_px:margin_px + length_px, x0:x0 + width_px] = 255
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(img).save(path)


def write_surface(path, gsd_m: float) -> None:
    """Minimal surface JSON in the form the pipeline reads: the stated scale (used as is) and
    a 1x1x3 float32 xyz.npy beside it."""
    path = Path(path)
    npy = path.with_name(path.stem + "_xyz.npy")
    np.save(npy, np.zeros((1, 1, 3), dtype=np.float32))
    meta = {"render_size_px": 1, "meters_per_pixel_x": gsd_m, "meters_per_pixel_y": gsd_m,
            "mapping_outputs": {"xyz_npy": npy.name}}
    path.write_text(json.dumps(meta, indent=2), encoding="utf-8")


def default_config_path() -> Path:
    """$CRACK_PIPELINE_CONFIG, else <repo root>/pipeline.local.json."""
    env = os.environ.get("CRACK_PIPELINE_CONFIG")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "pipeline.local.json"


def load_pipeline_config(path=None) -> dict:
    """{"cwd": directory to run in, "argv": the command up to --mask/--surface} (see
    pipeline.example.json)."""
    p = Path(path) if path is not None else default_config_path()
    if not p.is_file():
        raise FileNotFoundError(f"pipeline configuration not found at {p}; copy pipeline.example.json "
                                "to pipeline.local.json and fill it in")
    return json.loads(p.read_text(encoding="utf-8"))


def pipeline_command(config, mask, surface) -> list[str]:
    return [*config["argv"], "--mask", str(Path(mask).resolve()), "--surface", str(Path(surface).resolve())]


def run_pipeline(config, mask, surface) -> str:
    cp = subprocess.run(pipeline_command(config, mask, surface), cwd=config["cwd"],
                        capture_output=True, text=True, encoding="utf-8", check=True)
    return cp.stdout


_GSD = re.compile(r"GSD:\s+([\d.]+) mm/px\s+\(min measurable width ([\d.]+) mm\)")
_CRACK = re.compile(r"--- crack (\d+): (\d+) pts")
_WIDTH = re.compile(r"max width: ([\d.]+) px = ([\d.]+) ± ([\d.]+) mm")
_RATING = re.compile(r"DOH rating \(rule engine\): (.+?) — (.+)$")
_REFUSED = re.compile(r"REFUSED: (\S+)")


def parse_pipeline_output(text: str) -> dict:
    m = _GSD.search(text)
    if not m:
        raise ValueError("no GSD line in the pipeline output")
    out = {"gsd_mm": float(m[1]), "min_measurable_mm": float(m[2]), "cracks": []}
    cur = None
    for line in text.splitlines():
        if c := _CRACK.search(line):
            cur = {"index": int(c[1]), "points": int(c[2]), "gate": "ok", "max_width_px": None,
                   "max_width_mm": None, "uncertainty_mm": None, "rating": None, "rating_label": None}
            out["cracks"].append(cur)
        elif cur is not None and (w := _WIDTH.search(line)):
            cur.update(max_width_px=float(w[1]), max_width_mm=float(w[2]), uncertainty_mm=float(w[3]))
        elif cur is not None and (r := _RATING.search(line)):
            cur.update(rating=r[1].strip(), rating_label=r[2].strip())
        elif cur is not None and (f := _REFUSED.search(line)):
            cur["gate"] = f[1]
    return out


def governing_width(parsed: dict) -> dict | None:
    """The measured crack with the largest width, or None if every crack was refused."""
    ok = [c for c in parsed["cracks"] if c["gate"] == "ok"]
    return max(ok, key=lambda c: c["max_width_mm"]) if ok else None
```

- [ ] **Step 5: Run unit tests, then the live test**

Run: `python -m pytest tests/test_claimbridge.py -v`
Expected: 4 passed, 1 skipped

Run (with `pipeline.local.json` filled in): `python -m pytest tests/test_claimbridge.py -v`
Expected: 5 passed

- [ ] **Step 6: Run ClaimBridge's own crack-width test suite locally (spec §6: before relying on it)**

Run: the pipeline's own test command for its crack-width tests, in its own repository (not recorded here).
Expected: 0 failed. Record that the suite passed and the date for the progress log (Task 22).

- [ ] **Step 7: Commit**

```bash
git add src/crackfem/claimbridge.py tests/test_claimbridge.py tests/fixtures/spike_ok.txt tests/fixtures/spike_refused.txt
git commit -m "feat: file-only bridge to the pipeline (masks, surface, parser)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 19: Synthetic scenarios through the real pipeline

**Files:**
- Create: `scripts/make_synthetic_cases.py`, `scripts/run_claimbridge.py`
- Create (generated, committed): `data/claimbridge/cases.json`, `data/claimbridge/case*_mask.png`, `*_surface.json`, `*_surface_xyz.npy`, `*_spike.json`

- [ ] **Step 1: Write `scripts/make_synthetic_cases.py`**

```python
"""Synthetic scenarios: eta_true -> predicted crack width -> masks at several GSDs.

GSDs: 0.02 and 0.05 mm/px (close-range photo), 0.10 mm/px, 0.38 mm/px (lower end of the
thesis pipeline's field GSD). Widths are rounded to whole pixels, as a camera would.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem.claimbridge import write_mask, write_surface  # noqa: E402
from crackfem.girder1d import Girder1D  # noqa: E402
from crackfem.inverse import eta_decompression  # noqa: E402

GSDS_MM = (0.02, 0.05, 0.10, 0.38)
DROPS = (0.10, 0.20)          # eta_true = eta_decompression - drop


def main() -> None:
    g = Girder1D()
    eta_dec = eta_decompression(g.w_pred)
    out = ROOT / "data" / "claimbridge"
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    for k, drop in enumerate(DROPS, start=1):
        eta_true = round(eta_dec - drop, 2)
        w_true_mm = g.w_pred(eta_true) * 1e3
        for gsd in GSDS_MM:
            width_px = max(1, round(w_true_mm / gsd))
            name = f"case{k}_gsd{gsd:.2f}".replace(".", "p")
            write_mask(out / f"{name}_mask.png", width_px)
            write_surface(out / f"{name}_surface.json", gsd * 1e-3)
            cases.append({"name": name, "eta_true": eta_true, "w_true_mm": w_true_mm,
                          "gsd_mm": gsd, "width_px": width_px})
    (out / "cases.json").write_text(json.dumps({"eta_decompression": eta_dec, "cases": cases}, indent=2),
                                    encoding="utf-8")
    for c in cases:
        print(f"{c['name']}: eta_true={c['eta_true']:.2f} w_true={c['w_true_mm']:.3f} mm -> {c['width_px']} px")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write `scripts/run_claimbridge.py`**

```python
"""Run the thesis pipeline on every synthetic case and save the parsed report.

Uses the per-machine pipeline configuration (pipeline.local.json; see pipeline.example.json).
Machine-specific absolute paths are removed from the saved raw report before it is committed.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem.claimbridge import load_pipeline_config, parse_pipeline_output, run_pipeline  # noqa: E402


def main() -> None:
    config = load_pipeline_config()
    data = ROOT / "data" / "claimbridge"
    cases = json.loads((data / "cases.json").read_text(encoding="utf-8"))["cases"]
    for c in cases:
        text = run_pipeline(config, data / f"{c['name']}_mask.png", data / f"{c['name']}_surface.json")
        raw = "\n".join(l for l in text.splitlines() if not l.startswith(("Mask:", "Surface:")))
        record = {"case": c, "parsed": parse_pipeline_output(text), "raw": raw}
        (data / f"{c['name']}_spike.json").write_text(json.dumps(record, indent=2, ensure_ascii=False),
                                                     encoding="utf-8")
        print(c["name"], "->", [(k["gate"], k["max_width_mm"]) for k in record["parsed"]["cracks"]])


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Generate the masks**

Run: `python scripts/make_synthetic_cases.py`
Expected: 8 lines `caseK_gsd0pXX: eta_true=... w_true=... mm -> N px`.

- [ ] **Step 4: Run the real pipeline**

Run (with `pipeline.local.json` filled in): `python scripts/run_claimbridge.py`
Expected: 8 lines. Cases with fewer than 3 px show `('below-measurable-resolution', None)`; the others show `('ok', <width>)` for each of the 3 cracks.

- [ ] **Step 5: Check no private path leaked into committed data**

Run: `grep -rn "Users" data/claimbridge/*.json || echo "clean"`
Expected: `clean`

- [ ] **Step 6: Commit**

```bash
git add scripts/make_synthetic_cases.py scripts/run_claimbridge.py data/claimbridge
git commit -m "data: synthetic crack scenarios measured by the thesis pipeline" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 20: `figures.py`

**Files:**
- Create: `src/crackfem/figures.py`
- Test: `tests/test_figures.py`

- [ ] **Step 1: Write the failing test** (smoke: every function writes a non-empty PNG from small arrays)

```python
from crackfem import figures as F


def test_every_figure_function_writes_a_png(tmp_path):
    case = {"name": "c", "gsd_mm": 0.05, "eta_true": 0.4, "status": "measured", "w_meas_mm": 0.15,
            "u95_mm": 0.10, "eta_point": 0.41, "eta_mc": {"p05": 0.3, "p50": 0.41, "p95": 0.5}}
    refused = {**case, "name": "r", "gsd_mm": 0.38, "status": "refused", "eta_point": None, "eta_mc": None}
    calls = [
        (F.fig_w_eta, ([0.1, 0.5, 1.0], [0.4, 0.1, 0.0], [case])),
        (F.fig_mkappa, ([0, 1000, 2000], {"uncracked": [0, 1e-4, 2e-4], "cracked": [0, 2e-4, 5e-4]})),
        (F.fig_ei_ratio, ([1.0, 2.0, 3.0], {"eta = 1": [1, 1, 1], "eta = 0.4": [1, 0.6, 1]})),
        (F.fig_deflection, ([0, 1, 2], {"eta = 1": [0, 10, 0], "eta = 0.4": [0, 14, 0]})),
        (F.fig_eta_hist, ([0.3, 0.4, 0.41, 0.5], 0.4)),
        (F.fig_mesh, ([10, 20, 40], {"intact": [5, 5, 5], "damaged": [7, 7.2, 7.25]})),
        (F.fig_cases, ([case, refused],)),
        (F.fig_tornado, ([{"name": "a", "low": 1, "high": 2, "base": 0.4, "metric_low": 0.35, "metric_high": 0.45}],)),
    ]
    for i, (fn, args) in enumerate(calls):
        p = tmp_path / f"f{i}.png"
        fn(*args, p)
        assert p.stat().st_size > 1000, fn.__name__
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_figures.py -v`
Expected: FAIL — `cannot import name 'figures'`

- [ ] **Step 3: Write `src/crackfem/figures.py`**

```python
"""Report figures. Every function takes plain arrays/dicts (from results JSON) and a path."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _save(fig, path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_w_eta(eta, w_mm, cases, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(eta, w_mm, color="k", label="predicted crack width at rest")
    for c in cases:
        if c["status"] != "measured" or c["eta_point"] is None:
            continue
        ax.axhspan(c["w_meas_mm"] - c["u95_mm"], c["w_meas_mm"] + c["u95_mm"], alpha=0.12)
        ax.plot([c["eta_point"]], [c["w_meas_mm"]], "o", label=f"{c['name']} (GSD {c['gsd_mm']} mm/px)")
    ax.set_xlabel("remaining prestress fraction η")
    ax.set_ylabel("crack width w (mm)")
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_mkappa(M_kNm, curves, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, k in curves.items():
        ax.plot([1e3 * v for v in k], M_kNm, label=label)
    ax.set_xlabel("curvature from the end of stage A (1/km)")
    ax.set_ylabel("moment (kN·m)")
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_ei_ratio(x, ratios, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for label, r in ratios.items():
        ax.step(x, r, where="mid", label=label)
    ax.set_xlabel("x along the span (m)")
    ax.set_ylabel("secant EI / uncracked EI")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_deflection(x, shapes, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for label, v in shapes.items():
        ax.plot(x, v, label=label)
    ax.invert_yaxis()
    ax.set_xlabel("x along the span (m)")
    ax.set_ylabel("deflection, stage-B dead + truck (mm, down)")
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_eta_hist(samples, eta_true, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(samples, bins=40, color="0.6")
    ax.axvline(eta_true, color="k", linestyle="--", label="η true")
    ax.set_xlabel("η estimate (Monte Carlo on w)")
    ax.set_ylabel("count")
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_mesh(Ns, series, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for label, v in series.items():
        ax.plot(Ns, v, "o-", label=label)
    ax.set_xscale("log")
    ax.set_xlabel("number of elements N")
    ax.set_ylabel("midspan deflection, dead + truck (mm)")
    ax.legend(fontsize=7)
    _save(fig, path)


def fig_cases(cases, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for i, c in enumerate(cases):
        ax.plot([i], [c["eta_true"]], "kx")
        if c["status"] == "measured" and c["eta_mc"]:
            m = c["eta_mc"]
            ax.errorbar([i], [m["p50"]], yerr=[[m["p50"] - m["p05"]], [m["p95"] - m["p50"]]], fmt="o", capsize=3)
        else:
            ax.annotate("refused", (i, c["eta_true"]), textcoords="offset points", xytext=(0, 8),
                        ha="center", fontsize=7)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([f"{c['name']}" for c in cases], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("η (x = true, o = estimate 5–95%)")
    _save(fig, path)


def fig_tornado(rows, path) -> None:
    fig, ax = plt.subplots(figsize=(6, 0.5 + 0.45 * len(rows)))
    for i, r in enumerate(reversed(rows)):
        lo, hi = sorted((r["metric_low"], r["metric_high"]))
        ax.barh(i, hi - lo, left=lo, color="0.6")
    ax.axvline(rows[0]["base"], color="k")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["name"] for r in reversed(rows)], fontsize=8)
    ax.set_xlabel("η estimate")
    _save(fig, path)
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_figures.py -v`
Expected: 1 passed

- [ ] **Step 5: Commit**

```bash
git add src/crackfem/figures.py tests/test_figures.py
git commit -m "feat: report figures from plain arrays" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 21: `scripts/run_all.py` — results, convergence, cross-check, timing

**Files:**
- Create: `scripts/run_all.py`
- Create (generated, committed): `results/stage1.json`, `results/figures/*.png`

- [ ] **Step 1: Write `scripts/run_all.py`**

```python
"""Regenerate every Stage-1 result and figure from committed inputs.

Reads data/claimbridge/*_spike.json (thesis-pipeline output) and writes
results/stage1.json + results/figures/*.png. Does not need the thesis pipeline.
"""
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem import figures as F  # noqa: E402
from crackfem.claimbridge import governing_width  # noqa: E402
from crackfem.girder1d import Girder1D  # noqa: E402
from crackfem.inverse import estimate_eta, eta_decompression, is_monotone_nonincreasing  # noqa: E402
from crackfem.params import Inputs  # noqa: E402
from crackfem.solve import solve  # noqa: E402
from crackfem.uq import monte_carlo_eta, tornado  # noqa: E402

MC_N = 2000
MESH_NS = (10, 20, 40, 80, 160)
TORNADO = {"girder_spacing": (1.8, 2.6), "fck_girder": (35e6, 45e6), "strand_cover": (0.035, 0.055),
           "strand_centroid": (0.08, 0.12), "slab_thickness": (0.18, 0.22), "bond_xi": (0.5, 0.7)}


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def run_case(g: Girder1D, rec: dict) -> dict:
    c = rec["case"]
    gov = governing_width(rec["parsed"])
    out = {"name": c["name"], "gsd_mm": c["gsd_mm"], "eta_true": c["eta_true"], "w_true_mm": c["w_true_mm"],
           "width_px": c["width_px"]}
    if gov is None:
        # A refusal means w < 3 px, so eta lies above the eta whose predicted width is 3 px.
        gate = estimate_eta(g.w_pred, 3.0 * c["gsd_mm"] * 1e-3)
        return {**out, "status": "refused", "eta_point": None, "eta_mc": None,
                "reason": rec["parsed"]["cracks"][0]["gate"],
                "eta_lower_bound_from_gate": gate.eta if gate.status == "ok" else None,
                "gate_inverse_status": gate.status}
    w, u = gov["max_width_mm"] * 1e-3, gov["uncertainty_mm"] * 1e-3
    est = estimate_eta(g.w_pred, w)
    t0 = time.perf_counter()
    mc = monte_carlo_eta(g.w_pred, w, u, n=MC_N)
    t_mc = time.perf_counter() - t0
    s = mc.summary()
    band = {k: g.indicator(s[k]) for k in ("p05", "p50", "p95")} if s else {}
    return {**out, "status": "measured", "w_meas_mm": gov["max_width_mm"], "u95_mm": gov["uncertainty_mm"],
            "rating": gov["rating"], "rating_label": gov["rating_label"],
            "eta_point": est.eta, "eta_status": est.status, "eta_mc": s, "mc_status_counts": mc.status_counts,
            "mc_seconds": t_mc, "mc_samples": mc.etas,
            "deflection_service_mm_at": {k: v["deflection_service_mm"] for k, v in band.items()},
            "deflection_truck_mm_at": {k: v["deflection_truck_mm"] for k, v in band.items()}}


def main() -> None:
    res: dict = {}
    g = Girder1D()
    res["inputs"] = asdict(Inputs())
    res["design"] = asdict(g.design)
    res["girder"] = {"q_stage_a_kN_m": g.q_a / 1e3, "q_stage_b_kN_m": g.q_b / 1e3,
                     "M_A_mid_kNm": float(g.M_A_nd[g.mid_node]) / 1e3, "M_dead_mid_kNm": g.M_dead_mid / 1e3,
                     "M_service_kNm": g.M_service / 1e3, "EI_I_GNm2": g.EI_I / 1e9,
                     "EI_girder_GNm2": g.EI_girder / 1e9, "truck": g.truck, "h_m": g.h}
    grid = np.linspace(0.02, 1.0, 99)
    res["w_curve"] = {"eta": grid, "w_mm": [g.w_pred(float(e)) * 1e3 for e in grid],
                      "monotone": is_monotone_nonincreasing(g.w_pred, grid)}
    res["eta_decompression"] = eta_decompression(g.w_pred)

    data = ROOT / "data" / "claimbridge"
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(data.glob("*_spike.json"))]
    res["cases"] = [run_case(g, r) for r in recs]

    measured = [c for c in res["cases"] if c["status"] == "measured"]
    best = min(measured, key=lambda c: (c["gsd_mm"], -c["eta_true"]))
    eta_ref = best["eta_true"]            # finest GSD, least damage: the first detectable case
    res["eta_ref"] = eta_ref
    res["indicator_intact"] = g.indicator(1.0)
    res["indicator_ref"] = g.indicator(eta_ref)
    res["indicator_ref_branson"] = g.indicator(eta_ref, method="branson")
    M_grid = np.linspace(float(g.M_A_nd[g.mid_node]), 1.3 * g.M_service, 60)
    res["mkappa"] = {"M_kNm": M_grid / 1e3, **g.kappa_curves(eta_ref, M_grid)}

    res["mesh"] = {"N": list(MESH_NS),
                   "intact_mm": [Girder1D(n_el=n).indicator(1.0)["deflection_service_mm"] for n in MESH_NS],
                   "damaged_mm": [Girder1D(n_el=n).indicator(eta_ref)["deflection_service_mm"] for n in MESH_NS]}

    w_best = best["w_meas_mm"] * 1e-3
    res["mc_stability"] = {n: monte_carlo_eta(g.w_pred, w_best, best["u95_mm"] * 1e-3, n=n).summary()
                           for n in (1000, 10000)}

    def metric(inp: Inputs) -> float:
        est = estimate_eta(Girder1D(inp).w_pred, w_best)
        return est.eta if est.status == "ok" else float("nan")

    res["tornado"] = {"case": best["name"], "rows": tornado(metric, Inputs(), TORNADO)}

    t0 = time.perf_counter(); solve(g.model_service, g.EI_I); t_lin = time.perf_counter() - t0
    t0 = time.perf_counter(); g.analyse(eta_ref, g.model_service); t_nl = time.perf_counter() - t0
    t0 = time.perf_counter(); estimate_eta(g.w_pred, w_best); t_inv = time.perf_counter() - t0
    res["timing_s"] = {"linear_solve": t_lin, "nonlinear_analysis": t_nl, "inverse": t_inv,
                       "monte_carlo_best_case": best["mc_seconds"], "mc_samples": MC_N}

    out = ROOT / "results"
    figs = out / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    (out / "stage1.json").write_text(json.dumps(jsonable(res), indent=2, ensure_ascii=False), encoding="utf-8")

    F.fig_w_eta(res["w_curve"]["eta"], res["w_curve"]["w_mm"], measured, figs / "fig1_w_eta.png")
    F.fig_mkappa(res["mkappa"]["M_kNm"], {k: v for k, v in res["mkappa"].items() if k != "M_kNm"},
                 figs / "fig2_mkappa.png")
    F.fig_ei_ratio(res["indicator_ref"]["x_mid"], {"η = 1": res["indicator_intact"]["EI_ratio"],
                                                   f"η = {eta_ref}": res["indicator_ref"]["EI_ratio"]},
                   figs / "fig3_ei_ratio.png")
    F.fig_deflection(res["indicator_ref"]["x_nodes"], {"η = 1": res["indicator_intact"]["v_service_mm"],
                                                       f"η = {eta_ref}": res["indicator_ref"]["v_service_mm"]},
                     figs / "fig4_deflection.png")
    F.fig_eta_hist(best["mc_samples"], best["eta_true"], figs / "fig5_eta_hist.png")
    F.fig_mesh(res["mesh"]["N"], {"intact": res["mesh"]["intact_mm"], f"η = {eta_ref}": res["mesh"]["damaged_mm"]},
               figs / "fig6_mesh.png")
    F.fig_cases(res["cases"], figs / "fig7_cases.png")
    F.fig_tornado(res["tornado"]["rows"], figs / "fig8_tornado.png")

    print(f"eta_decompression = {res['eta_decompression']:.4f}")
    for c in res["cases"]:
        print(c["name"], c["status"], c.get("eta_point"), c.get("eta_mc"))
    print("timing", res["timing_s"])


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

Run: `python scripts/run_all.py`
Expected: prints `eta_decompression = ...`, one line per case (refused cases at the coarse GSDs, measured ones with an η estimate near `eta_true`), and the timing dict. Exit code 0.

- [ ] **Step 3: Check the results against the tests' expectations (not by eye only)**

Run: `python -m pytest -q`
Expected: all passed. Then open `results/stage1.json` and confirm, writing each into the log in Task 22:
- `w_curve.monotone` is `true`
- `indicator_ref.iterations_service` is `1`
- `mesh.intact_mm` values are equal to each other (nodal exactness) and `mesh.damaged_mm` changes by less than 1% from N = 80 to 160
- every `cases[*].mc_status_counts` has only `"ok"` for measured cases

If any of these fails, stop and debug (superpowers:systematic-debugging) before Task 22.

- [ ] **Step 4: Look at every figure** (`results/figures/fig1…fig8.png`) with the Read tool, one by one; fix labels or scales in `figures.py` if a figure is unreadable, then re-run Step 2.

- [ ] **Step 5: Commit**

```bash
git add scripts/run_all.py results
git commit -m "results: Stage-1 run, mesh convergence, Branson cross-check, timing, figures" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Task 22: README, progress log B3, reading guide, final verification

**Files:**
- Create: `README.md`
- Modify: `docs/PROGRESS_LOG.md`, `docs/DECISIONS.md`

- [ ] **Step 1: Write `README.md`**

```markdown
# crackfem — crack width → remaining prestress → stiffness (FEM course project)

01203541 Finite Element Analysis of Structures, Kasetsart University, 2026 · Palid Archavachalee (6914500959)

A measured crack width on a prestressed bridge girder is turned into a mechanical statement:
how much prestress is left (η), how stiff the girder still is along the span, and how far it
deflects under the design truck — with an uncertainty band from the measurement. Crack widths
come from the author's thesis inspection pipeline (private) run on synthetic crack images.
Stage 1 (this release) models one girder; Stage 2 adds a grillage of the whole deck.

## Quickstart

    python -m pip install -r requirements.txt
    python -m pytest -q
    python scripts/run_all.py

`run_all.py` regenerates `results/stage1.json` and `results/figures/` from the committed
pipeline outputs in `data/claimbridge/`; it does not need the private repository.

## Where to read

- Design and the five changes from the proposal: `docs/superpowers/specs/`
- What happened, block by block: `docs/PROGRESS_LOG.md`
- Every input value and its source: `docs/sources.md`
- Decisions and rejected alternatives: `docs/DECISIONS.md`

## Authorship

The code was written by Claude (Anthropic's AI) and reviewed by the student; review notes
are in `docs/PROGRESS_LOG.md`.

## Limits

All structural values are stand-ins until the as-built drawings are available, and crack
widths are synthetic until the November 2026 field survey (see the spec §5.6).
```

- [ ] **Step 2: Verify the Ch3 / Ch11 references for the reading guide.** Run `pdftotext` on the local copy of Logan (5th ed.) and grep for the beam stiffness matrix section and work-equivalence section; record book page numbers (PDF page = book page + 19). Use only page numbers found this way in Step 3.

- [ ] **Step 3: Append the B3 entry to `docs/PROGRESS_LOG.md`** (every `<...>` filled from `results/stage1.json`, test output and Task 18 Step 6; nothing in angle brackets may remain):

```markdown
## Entry 3 — B3 thesis-pipeline link, Monte Carlo, Stage-1 results (<date>)

**Planned:** run synthetic cracks through the real thesis pipeline, invert for η with a Monte Carlo band, verify mesh convergence and the Branson cross-check (spec §5.4–§6).

**Found:**
- Design: <n_strands> strands, P0 = <P0> kN, decompression η = <eta_decompression>.
- The thesis pipeline refuses cracks narrower than 3 px: at GSD <list> mm/px every case was "<reason>" — at field GSD it cannot see the crack widths that matter for a prestressed girder.
- Measured cases: <per case: GSD, w ± u, η estimate 5–50–95% vs η true, DOH rule-engine rating>.
- The rule engine rates <width> mm as "<rating label>", while the model says η ≈ <p50> (<100·(1−p50)>% of prestress lost).
- Mesh: intact deflection identical for N = 10…160 (nodal exactness); damaged deflection changes <x>% from N = 80 to 160.
- Visible crack zone at rest at η = <eta_ref>: <z> m. On synthetic data this is a model output only (the masks are drawn from the same model); it becomes the spec's held-out check once field images exist.
- Secant loop: 1 pass on the determinate girder (P5 confirmed) — the grillage (Plan 2) is where it matters.
- Branson vs EC2 service deflection at η = <eta_ref>: <a> vs <b> mm.
- Timing: linear solve <t> ms, nonlinear analysis <t> ms, inverse <t> ms, Monte Carlo (2,000) <t> s → no surrogate needed.
- ClaimBridge's own crack-width test suite: <passed, run locally on <date>>.

**Changed & why:** <deviations from Plan 1, if any>.

**Numbers:** `results/stage1.json`; figures `results/figures/fig1…fig8.png`; `python -m pytest -q` → <count> passed, <count> skipped.

**Reading guide:**

| Code | Course | Where |
|---|---|---|
| `elements.beam_k`, `hermite` | beam element stiffness from Hermite shape functions | Logan p. <page> |
| `elements.udl_equivalent`, `point_equivalent` | work-equivalent nodal loads (`f_eq = −F0` [superseded: f_eq = +{f0}, see As-built notes] of the slides) | Logan p. <page> |
| `solve.solve` | assembly + partition K_ff D_f = F_f, reactions | Ch3 slides / Logan p. <page> |
| `nonlinear.solve_secant` | direct iteration with secant stiffness | Ch11 slides |
| `rootfind.illinois` | secant method kept inside a bracket | Ch11 (secant idea) |

Explain-back questions: (1) Why does a 0.15 mm crack at rest mean prestress was lost, when the same crack in an RC beam is harmless? (2) Why can the field GSD not support this indicator? (3) Why is the Monte Carlo band wide, and what single change to the measurement would narrow it most?

**Reviewed by student:** pending (after 2026-10-19).
```

- [ ] **Step 4: Append to `docs/DECISIONS.md`** one line per non-obvious choice made during Plan 1 (at least: pre-cracked zone rule — elements whose service moment exceeds Mcr(η) lose tensile strength; dead-load midspan moment computed once on the determinate girder; Illinois as the safeguarded secant; GSD set 0.02/0.05/0.10/0.38 mm/px; held-out cracked-zone check deferred to field data because synthetic masks come from the same model; cracked-section check uses the textbook RC closed form plus the S21 prestressed example if one was found).

- [ ] **Step 5: Final verification**

Run: `python -m pytest -q` and `python -m pytest tests/test_claimbridge.py -q`
Expected: all passed (only the S21 skip allowed in the first run); 5 passed in the second. Paste both summary lines into the progress-log entry.

Run: `git status --short`
Expected: only the files of this task are modified.

- [ ] **Step 6: Commit**

```bash
git add README.md docs/PROGRESS_LOG.md docs/DECISIONS.md
git commit -m "docs: README, Stage-1 progress entry with reading guide" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

**Stage 1 is the fallback deliverable from here on.** Do not push to GitHub (gate: ask the student, including public vs private).

---

## Plan 2 scope (written 2026-10-19/20, not part of this plan)

- **B4 (Oct 20–23):** grid element (torsion + bending + [T]) reproducing a Logan §5.4 example; deck grillage model of S01 (girders + slab strips + diaphragms per Hambly rules from S17/sources); PyNite install (approved) and cross-check; grillage nonlinear loop (with under-relaxation if direct iteration oscillates).
- **B5 (Oct 24–26):** w_pred(η) through the nonlinear grillage (cracked girder G1, second case G3); load shed per girder; Monte Carlo and tornado on the grillage; final figures.
- **B6 (Oct 27–30):** report (docx + PDF, English) assembled from `PROGRESS_LOG.md`; student review sessions; fixes; GitHub push after the gate.
