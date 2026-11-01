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
import shutil
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image


def write_mask(path, width_px: int, n_cracks: int = 3, length_px: int = 400,
               spacing_px: int = 150, margin_px: int = 50) -> None:
    """Vertical flexural cracks of constant width, white on black."""
    if isinstance(width_px, bool) or not isinstance(width_px, (int, np.integer)) or width_px < 1:
        raise ValueError("width_px must be an integer of at least 1")
    if n_cracks < 1:
        raise ValueError("n_cracks must be at least 1")
    if width_px >= spacing_px:
        raise ValueError("width_px must be smaller than spacing_px, or neighbouring cracks merge")
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
    if not (np.isfinite(gsd_m) and gsd_m > 0):
        raise ValueError("gsd_m must be a positive, finite number (metres per pixel)")
    path = Path(path)
    npy = path.with_name(path.stem + "_xyz.npy")
    np.save(npy, np.zeros((1, 1, 3), dtype=np.float32))
    meta = {"render_size_px": 1, "meters_per_pixel_x": gsd_m, "meters_per_pixel_y": gsd_m,
            "mapping_outputs": {"xyz_npy": npy.name}}
    path.write_text(json.dumps(meta, indent=2), encoding="utf-8", newline="\n")


CONFIG_ENV = "CRACK_PIPELINE_CONFIG"
CONFIG_NAME = "pipeline.local.json"
EXAMPLE_NAME = "pipeline.example.json"
_REPO_ROOT = Path(__file__).resolve().parents[2]


def default_config_path() -> Path:
    """Where the pipeline configuration is looked up when no path is given: the file named by
    the CRACK_PIPELINE_CONFIG environment variable, else pipeline.local.json in the repo root."""
    env = os.environ.get(CONFIG_ENV)
    return Path(env) if env else _REPO_ROOT / CONFIG_NAME


def load_pipeline_config(path=None) -> dict:
    """Read the per-machine pipeline configuration: {"cwd": directory to run in, "argv": the
    command up to, not including, the --mask/--surface options}.

    `path` if given, else CRACK_PIPELINE_CONFIG, else <repo root>/pipeline.local.json.
    Raises FileNotFoundError (naming pipeline.example.json) when the file is missing and
    ValueError when it is not valid JSON or lacks a usable "cwd" / "argv"."""
    p = Path(path) if path is not None else default_config_path()
    if not p.is_file():
        raise FileNotFoundError(
            f"pipeline configuration not found at {p}; copy {EXAMPLE_NAME} to {CONFIG_NAME} "
            f"and fill it in (or point {CONFIG_ENV} to your file)")
    try:
        cfg = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise ValueError(f"pipeline configuration {p} is not readable JSON: {e}") from e
    if not isinstance(cfg, dict):
        raise ValueError(f"pipeline configuration {p} must be a JSON object")
    cwd, argv = cfg.get("cwd"), cfg.get("argv")
    if not isinstance(cwd, str) or not cwd.strip():
        raise ValueError(f'pipeline configuration {p}: "cwd" must be a non-empty string')
    if (not isinstance(argv, list) or not argv
            or not all(isinstance(a, str) and a.strip() for a in argv)):
        raise ValueError(f'pipeline configuration {p}: "argv" must be a non-empty list of strings')
    return {"cwd": cwd, "argv": list(argv)}


def pipeline_command(config: dict, mask, surface) -> list[str]:
    """The full command line: the configured argv, then --mask and --surface as absolute paths
    (resolved here, against this process's working directory, because the pipeline runs in
    its own configured directory)."""
    return [*config["argv"], "--mask", str(Path(mask).resolve()), "--surface", str(Path(surface).resolve())]


def run_pipeline(config: dict, mask, surface, timeout: float = 120) -> str:
    """Run the configured pipeline command on a mask and a surface file and return its stdout.

    Raises RuntimeError (carrying the pipeline's stderr) when the executable is missing, the
    command exits non-zero, or it runs longer than `timeout` seconds; FileNotFoundError when
    the configured directory does not exist."""
    argv = config["argv"]
    if shutil.which(argv[0]) is None:
        raise RuntimeError(f"{argv[0]!r} was not found on PATH; check argv in {CONFIG_NAME}")
    cwd = Path(config["cwd"])
    if not cwd.is_dir():
        raise FileNotFoundError(f"pipeline directory not found at {cwd}; check cwd in {CONFIG_NAME}")
    cmd = pipeline_command(config, mask, surface)
    try:
        cp = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, encoding="utf-8",
                            check=True, timeout=timeout)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"the pipeline failed (exit code {e.returncode}): {e.stderr}") from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(f"the pipeline timed out after {timeout} s") from e
    return cp.stdout


_GSD = re.compile(r"GSD:\s+([\d.]+) mm/px\s+\(min measurable width ([\d.]+) mm\)")
_CRACK_START = re.compile(r"^\s*--- crack\b")
_CRACK = re.compile(r"--- crack (\d+): (\d+) pts(?: \((\d+) refused, (\d+) flagged\))?")
_WIDTH = re.compile(r"max width: ([\d.]+) px = ([\d.]+) ± ([\d.]+) mm")
_RATING = re.compile(r"DOH rating \(rule engine\): (.+?) — (.+)$")
_REFUSED = re.compile(r"REFUSED: (\S+)")

# A crack refused with one of these gates was read, but it is narrower than the minimum
# measurable width, so it cannot be wider than a measured crack. Any other non-"ok" gate
# (notably "no-valid-points") leaves the crack's width unknown.
_NARROWER_THAN_MEASURABLE = ("below-measurable-resolution",)


def recognised_lines(text: str) -> list[str]:
    """The lines of the pipeline's report that parse_pipeline_output understands (this drops
    the Mask:/Surface: lines, which carry local file paths)."""
    return [ln for ln in text.splitlines()
            if any(rx.search(ln) for rx in (_GSD, _CRACK, _WIDTH, _RATING, _REFUSED))]


def parse_pipeline_output(text: str) -> dict:
    """Parse the pipeline's report. Fails closed: a crack whose result is neither a width line
    nor a REFUSED line (e.g. a changed output format) raises ValueError instead of passing."""
    m = _GSD.search(text)
    if not m:
        raise ValueError("no GSD line in the pipeline output")
    out = {"gsd_mm": float(m[1]), "min_measurable_mm": float(m[2]), "cracks": []}
    cur = None
    for line in text.splitlines():
        if _CRACK_START.search(line):
            c = _CRACK.search(line)
            if not c:
                raise ValueError(f"unrecognised crack header in the pipeline output: {line.strip()!r}")
            cur = {"index": int(c[1]), "points": int(c[2]),
                   "refused_points": None if c[3] is None else int(c[3]),
                   "flagged_points": None if c[4] is None else int(c[4]),
                   "gate": None, "max_width_px": None, "max_width_mm": None, "uncertainty_mm": None,
                   "rating": None, "rating_label": None}
            out["cracks"].append(cur)
        elif cur is not None and (w := _WIDTH.search(line)):
            cur.update(gate="ok", max_width_px=float(w[1]), max_width_mm=float(w[2]),
                       uncertainty_mm=float(w[3]))
        elif cur is not None and (r := _RATING.search(line)):
            cur.update(rating=r[1].strip(), rating_label=r[2].strip())
        elif cur is not None and (f := _REFUSED.search(line)):
            cur["gate"] = f[1]
    for c in out["cracks"]:
        if c["gate"] is None:
            raise ValueError(f"crack {c['index']}: neither a width line nor a REFUSED line in the output")
    return out


def governing_width(parsed: dict) -> dict | None:
    """The measured crack with the largest width, or None if there is no safe governing width.

    None when every crack was refused, and also when any crack is unresolved (gate
    "no-valid-points", or a gate this module does not know): an unresolved crack may be the
    widest one, so reporting the widest *measured* crack would bias the width low, i.e. the
    remaining prestress eta high, which is the unconservative direction. A crack refused only
    for being narrower than the minimum measurable width cannot exceed the measured ones and
    is ignored."""
    cracks = parsed["cracks"]
    if any(c["gate"] not in ("ok", *_NARROWER_THAN_MEASURABLE) for c in cracks):
        return None
    ok = [c for c in cracks if c["gate"] == "ok"]
    return max(ok, key=lambda c: c["max_width_mm"]) if ok else None


def refused_below_resolution(parsed: dict) -> bool:
    """True only if there is at least one crack and every crack was refused as narrower than
    the minimum measurable width.

    Only then does the refusal bound the governing width from above (below the gate, i.e. 3 px),
    which is what an eta lower bound needs. A report with a measured ("ok") crack, with an
    unresolved one ("no-valid-points" or a gate this module does not know), or with no cracks at
    all gives no such bound: an unresolved crack may be wider than the gate."""
    cracks = parsed["cracks"]
    return bool(cracks) and all(c["gate"] in _NARROWER_THAN_MEASURABLE for c in cracks)
