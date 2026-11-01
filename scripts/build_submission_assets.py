"""Build compact submission assets from the committed Stage-1 and Stage-2 results.

Writes a workflow figure, a flat key-results table, and a provenance record. It does not
recompute the FEM models or modify stage1.json/stage2.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def git(*args: str) -> str:
    cp = subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True)
    return cp.stdout.strip()


def workflow_figure(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 3.1), constrained_layout=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.02, "Crack image\nand surface scale", "Field / synthetic input"),
        (0.215, "ClaimBridge", "Mask + GSD to width w +/- u"),
        (0.41, "Inverse problem", "Width to remaining prestress eta"),
        (0.605, "Nonlinear FEM", "Element EI and deck load sharing"),
        (0.80, "Engineering outputs", "Deflection, stiffness and moment"),
    ]
    width = 0.17
    for x, title, subtitle in boxes:
        patch = FancyBboxPatch((x, 0.27), width, 0.48, boxstyle="round,pad=0.012,rounding_size=0.02",
                               facecolor="#EAF2F8", edgecolor="#245A7A", linewidth=1.5)
        ax.add_patch(patch)
        ax.text(x + width / 2, 0.60, title, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x + width / 2, 0.40, subtitle, ha="center", va="center", fontsize=8.5, wrap=True)
    for i in range(len(boxes) - 1):
        x1 = boxes[i][0] + width
        x2 = boxes[i + 1][0]
        ax.add_patch(FancyArrowPatch((x1 + 0.006, 0.51), (x2 - 0.006, 0.51),
                                     arrowstyle="-|>", mutation_scale=13, color="#444444", linewidth=1.2))
    ax.text(0.50, 0.11,
            "Course-repository boundary: the private image-to-mask AI is outside this model; the tested interface starts at the mask and GSD.",
            ha="center", va="center", fontsize=9, color="#333333")
    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)


def key_results(stage1: dict, stage2: dict, path: Path) -> None:
    rows = []
    model_cases = (("Stage 1 single girder", stage1["cases"]),
                   ("Stage 2 G1 exterior", stage2["G1"]["cases"]),
                   ("Stage 2 G3 interior", stage2["G3"]["cases"]))
    for model, cases in model_cases:
        for case in cases:
            if case["status"] != "measured":
                continue
            mc = case["eta_mc"]
            rows.append({
                "model": model,
                "case": case["name"],
                "gsd_mm_per_px": case["gsd_mm"],
                "width_mm": case["w_meas_mm"],
                "uncertainty_plus_minus_mm": case["u95_mm"],
                "eta_point": case["eta_point"],
                "eta_p05": mc["p05"],
                "eta_p50": mc["p50"],
                "eta_p95": mc["p95"],
                "mc_samples_used": mc["n"],
                "mc_samples_dropped": mc["n_dropped"],
                "evidence_scope": "synthetic mask; measurement noise only; stand-in structure",
            })
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_record(stage1: dict, stage2: dict, path: Path) -> None:
    def sha256(file: Path) -> str:
        return hashlib.sha256(file.read_bytes()).hexdigest()

    status = git("status", "--short").splitlines()
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository_head": git("rev-parse", "HEAD"),
        "repository_head_short": git("rev-parse", "--short", "HEAD"),
        "working_tree_clean_at_record_time": not bool(status),
        "working_tree_entries_at_record_time": status,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "source_results": {
            "stage1": "results/stage1.json",
            "stage2": "results/stage2.json",
            "proposal_checks": "results/proposal_checks.json",
            "stage1_sha256": sha256(RESULTS / "stage1.json"),
            "stage2_sha256": sha256(RESULTS / "stage2.json"),
            "proposal_checks_sha256": sha256(RESULTS / "proposal_checks.json"),
            "stage1_eta_decompression": stage1["eta_decompression"],
            "stage1_eta_reference": stage1["eta_ref"],
            "stage2_model": stage2["model"],
            "stage2_strands_G1_to_G5": stage2["design"]["n_strands"],
        },
        "verification_captured_2026_11_01": [
            {"command": "python -m pytest -q", "result": "425 passed", "elapsed_s": 84.95,
             "private_pipeline_config_present": True},
            {"command": "python -m pytest -q (fresh clone of the published commit)", "result": "422 passed, 3 skipped",
             "elapsed_s": 87.01, "private_pipeline_config_present": False},
        ],
        "scope_warning": "Synthetic masks and stand-in structure; not field validation or a capacity assessment.",
    }
    path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    stage1 = json.loads((RESULTS / "stage1.json").read_text(encoding="utf-8"))
    stage2 = json.loads((RESULTS / "stage2.json").read_text(encoding="utf-8"))
    workflow_figure(RESULTS / "figures" / "fig0_workflow.png")
    key_results(stage1, stage2, RESULTS / "key_results.csv")
    run_record(stage1, stage2, RESULTS / "run_record.json")
    print("wrote fig0_workflow.png, key_results.csv, run_record.json")


if __name__ == "__main__":
    main()
