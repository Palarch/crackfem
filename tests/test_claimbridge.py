import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from crackfem.claimbridge import (CONFIG_ENV, default_config_path, governing_width, load_pipeline_config,
                                  parse_pipeline_output, pipeline_command, refused_below_resolution,
                                  run_pipeline, write_mask, write_surface)

ROOT = Path(__file__).resolve().parents[1]
FIX = Path(__file__).parent / "fixtures"

# The live tests need the per-machine pipeline configuration (CRACK_PIPELINE_CONFIG or
# pipeline.local.json); without it they are skipped. Everything else runs offline.
needs_pipeline = pytest.mark.skipif(
    not default_config_path().is_file(),
    reason="needs a pipeline config (CRACK_PIPELINE_CONFIG or pipeline.local.json)")

_spec = importlib.util.spec_from_file_location("run_claimbridge", ROOT / "scripts" / "run_claimbridge.py")
run_claimbridge = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_claimbridge)


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


@needs_pipeline
def test_live_pipeline_measures_a_six_pixel_crack(tmp_path):
    write_mask(tmp_path / "m.png", width_px=6)
    write_surface(tmp_path / "s.json", 5e-5)
    out = parse_pipeline_output(run_pipeline(load_pipeline_config(), tmp_path / "m.png", tmp_path / "s.json"))
    assert governing_width(out)["max_width_mm"] == pytest.approx(0.30)


# ---- fail-closed parsing -------------------------------------------------------------------

HDR = "Mask:      m.png\nSurface:   s.json\nGSD:       0.0500 mm/px  (min measurable width 0.15 mm)\n"
OK_CRACK = ("--- crack {i}: 394 pts (0 refused, 0 flagged)\n    max width: {px} px = {mm} ± 0.10 mm\n"
            "    DOH rating (rule engine): 5 — ดีมาก (Excellent)\n")
REFUSED_CRACK = "--- crack {i}: 390 pts (0 refused, 0 flagged)\n    REFUSED: {gate}\n"


def _records():
    return sorted((Path(__file__).resolve().parents[1] / "data" / "claimbridge").glob("*_spike.json"))


@pytest.mark.parametrize("path", _records(), ids=lambda p: p.stem)
def test_committed_raw_reports_reparse_to_the_saved_values(path):
    """The saved raw text (paths removed) is enough to reproduce every parsed field it records."""
    rec = json.loads(path.read_text(encoding="utf-8"))
    again = parse_pipeline_output(rec["raw"])
    parsed = rec["parsed"]
    assert again["gsd_mm"] == parsed["gsd_mm"] and again["min_measurable_mm"] == parsed["min_measurable_mm"]
    assert len(again["cracks"]) == len(parsed["cracks"]) > 0
    for new, saved in zip(again["cracks"], parsed["cracks"]):
        assert {k: new[k] for k in saved} == saved


def test_crlf_output_parses_like_lf_output():
    lf = (FIX / "spike_ok.txt").read_text(encoding="utf-8").replace("\r\n", "\n")
    assert parse_pipeline_output(lf.replace("\n", "\r\n")) == parse_pipeline_output(lf)


def test_refused_and_flagged_point_counts_are_kept():
    text = HDR + "--- crack 1: 400 pts (350 refused, 40 flagged)\n    max width: 6.00 px = 0.30 ± 0.10 mm\n"
    c = parse_pipeline_output(text)["cracks"][0]
    assert (c["points"], c["refused_points"], c["flagged_points"]) == (400, 350, 40)
    assert c["gate"] == "ok" and c["rating"] is None


def test_unknown_result_line_raises_instead_of_passing_silently():
    text = HDR + OK_CRACK.format(i=1, px="3.00", mm="0.15") \
        + "--- crack 2: 10 pts (0 refused, 0 flagged)\n    UNMEASURABLE: edge-clipped\n"
    with pytest.raises(ValueError, match="crack 2"):
        parse_pipeline_output(text)


def test_width_line_without_the_uncertainty_raises():
    text = HDR + "--- crack 1: 10 pts (0 refused, 0 flagged)\n    max width: 3.00 px = 0.15 mm\n"
    with pytest.raises(ValueError, match="crack 1"):
        parse_pipeline_output(text)


def test_changed_crack_header_raises():
    with pytest.raises(ValueError, match="crack header"):
        parse_pipeline_output(HDR + "--- crack 1: 10 points\n    REFUSED: below-measurable-resolution\n")


def test_cli_error_text_is_not_mistaken_for_a_report():
    with pytest.raises(ValueError, match="GSD"):
        parse_pipeline_output("error: cannot determine metric scale\n")


def test_no_valid_points_crack_blocks_the_governing_width():
    # The unresolved crack may be the widest: dropping it would bias the width low (eta high).
    text = HDR + OK_CRACK.format(i=1, px="3.00", mm="0.15") + REFUSED_CRACK.format(i=2, gate="no-valid-points")
    out = parse_pipeline_output(text)
    assert out["cracks"][1]["gate"] == "no-valid-points"
    assert governing_width(out) is None


def test_unknown_gate_also_blocks_the_governing_width():
    text = HDR + OK_CRACK.format(i=1, px="3.00", mm="0.15") + REFUSED_CRACK.format(i=2, gate="some-new-gate")
    assert governing_width(parse_pipeline_output(text)) is None


def test_crack_narrower_than_measurable_does_not_block_the_governing_width():
    text = (HDR + OK_CRACK.format(i=1, px="3.00", mm="0.15") + OK_CRACK.format(i=2, px="5.00", mm="0.25")
            + REFUSED_CRACK.format(i=3, gate="below-measurable-resolution"))
    assert governing_width(parse_pipeline_output(text))["max_width_mm"] == 0.25


def test_report_without_cracks_has_no_governing_width():
    assert governing_width(parse_pipeline_output(HDR)) is None


def test_refused_below_resolution_only_when_every_crack_is_below_the_gate():
    below = REFUSED_CRACK.format(i=1, gate="below-measurable-resolution")
    assert refused_below_resolution(parse_pipeline_output(
        HDR + below + REFUSED_CRACK.format(i=2, gate="below-measurable-resolution"))) is True


def test_refused_below_resolution_is_false_with_a_no_valid_points_crack():
    # an unresolved crack may be wider than the gate: no bound from the refusal
    text = (HDR + REFUSED_CRACK.format(i=1, gate="below-measurable-resolution")
            + REFUSED_CRACK.format(i=2, gate="no-valid-points"))
    assert refused_below_resolution(parse_pipeline_output(text)) is False


def test_refused_below_resolution_is_false_with_a_measured_crack():
    text = (HDR + REFUSED_CRACK.format(i=1, gate="below-measurable-resolution")
            + OK_CRACK.format(i=2, px="3.00", mm="0.15"))
    assert refused_below_resolution(parse_pipeline_output(text)) is False


def test_refused_below_resolution_is_false_with_an_unknown_gate():
    text = HDR + REFUSED_CRACK.format(i=1, gate="some-new-gate")
    assert refused_below_resolution(parse_pipeline_output(text)) is False


def test_refused_below_resolution_is_false_without_cracks():
    assert refused_below_resolution(parse_pipeline_output(HDR)) is False


# ---- input validation ----------------------------------------------------------------------

@pytest.mark.parametrize("kwargs", [
    {"width_px": 0}, {"width_px": -3}, {"width_px": 2.5}, {"width_px": 3.0}, {"width_px": True},
    {"width_px": 5, "n_cracks": 0},
    {"width_px": 150}, {"width_px": 200}, {"width_px": 40, "spacing_px": 40},
])
def test_write_mask_rejects_unusable_arguments(tmp_path, kwargs):
    with pytest.raises(ValueError):
        write_mask(tmp_path / "m.png", **kwargs)
    assert not (tmp_path / "m.png").exists()


@pytest.mark.parametrize("gsd", [0, -5e-5, float("nan"), float("inf")])
def test_write_surface_rejects_a_non_positive_gsd(tmp_path, gsd):
    with pytest.raises(ValueError):
        write_surface(tmp_path / "s.json", gsd)
    assert not (tmp_path / "s.json").exists()


def test_surface_json_is_written_with_unix_newlines(tmp_path):
    write_surface(tmp_path / "s.json", 5e-5)
    assert b"\r" not in (tmp_path / "s.json").read_bytes()


# ---- pipeline configuration ----------------------------------------------------------------

def _write_config(path, cwd="somewhere", argv=("tool", "run")):
    path.write_text(json.dumps({"cwd": cwd, "argv": list(argv)}), encoding="utf-8")
    return path


def test_load_pipeline_config_reads_cwd_and_argv(tmp_path):
    f = _write_config(tmp_path / "c.json", cwd=str(tmp_path), argv=["tool", "a", "b"])
    assert load_pipeline_config(f) == {"cwd": str(tmp_path), "argv": ["tool", "a", "b"]}


def test_config_path_precedence_argument_then_env_then_repo_file(tmp_path, monkeypatch):
    explicit = _write_config(tmp_path / "explicit.json", cwd="from-argument")
    via_env = _write_config(tmp_path / "env.json", cwd="from-env")
    root = tmp_path / "root"
    root.mkdir()
    _write_config(root / "pipeline.local.json", cwd="from-repo-file")
    monkeypatch.setattr("crackfem.claimbridge._REPO_ROOT", root)
    monkeypatch.setenv(CONFIG_ENV, str(via_env))
    assert load_pipeline_config(explicit)["cwd"] == "from-argument"
    assert load_pipeline_config()["cwd"] == "from-env"
    monkeypatch.delenv(CONFIG_ENV)
    assert load_pipeline_config()["cwd"] == "from-repo-file"


def test_missing_config_raises_file_not_found_naming_the_example(tmp_path, monkeypatch):
    monkeypatch.setattr("crackfem.claimbridge._REPO_ROOT", tmp_path)
    monkeypatch.delenv(CONFIG_ENV, raising=False)
    with pytest.raises(FileNotFoundError, match=r"pipeline\.example\.json"):
        load_pipeline_config()
    with pytest.raises(FileNotFoundError, match=r"pipeline\.example\.json"):
        load_pipeline_config(tmp_path / "nope.json")


def test_config_named_by_the_environment_but_missing_also_names_the_example(tmp_path, monkeypatch):
    monkeypatch.setenv(CONFIG_ENV, str(tmp_path / "nope.json"))
    with pytest.raises(FileNotFoundError, match=r"pipeline\.example\.json"):
        load_pipeline_config()


@pytest.mark.parametrize("text", [
    "not json at all", "[]", '"a string"', "{}",
    '{"cwd": "x"}', '{"argv": ["a"]}',
    '{"cwd": "", "argv": ["a"]}', '{"cwd": 3, "argv": ["a"]}',
    '{"cwd": "x", "argv": []}', '{"cwd": "x", "argv": "tool run"}',
    '{"cwd": "x", "argv": ["a", 3]}', '{"cwd": "x", "argv": ["a", ""]}',
])
def test_malformed_config_raises_value_error(tmp_path, text):
    f = tmp_path / "c.json"
    f.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        load_pipeline_config(f)


def test_committed_example_config_is_well_formed_and_holds_only_placeholders():
    cfg = load_pipeline_config(ROOT / "pipeline.example.json")
    placeholders = [cfg["cwd"], *cfg["argv"]]
    assert all(v.startswith("<") and v.endswith(">") for v in placeholders)
    assert not any(c in v for v in placeholders for c in ":\\")


def test_local_config_is_ignored_by_git():
    assert "pipeline.local.json" in (ROOT / ".gitignore").read_text(encoding="utf-8").split()


# ---- running the pipeline (offline: a stand-in command written by the test) -------------------

FAKE_PIPELINE = '''\
import json, os, sys, time
mode = sys.argv[1]
if mode == "report":
    print(json.dumps({"cwd": os.getcwd(), "args": sys.argv[2:]}))
elif mode == "fail":
    sys.stderr.write("boom: cannot read the surface\\n")
    sys.exit(3)
elif mode == "sleep":
    time.sleep(30)
'''


@pytest.fixture
def fake(tmp_path):
    """Pipeline configs whose command is a small Python script in its own directory."""
    d = tmp_path / "pipeline_dir"
    d.mkdir()
    (d / "fake_pipeline.py").write_text(FAKE_PIPELINE, encoding="utf-8")

    def make(mode):
        return {"cwd": str(d), "argv": [sys.executable, "fake_pipeline.py", mode]}
    make.dir = d
    return make


def test_pipeline_command_appends_mask_and_surface_after_the_configured_argv(tmp_path):
    cfg = {"cwd": "unused", "argv": ["tool", "sub", "cmd"]}
    cmd = pipeline_command(cfg, tmp_path / "m.png", tmp_path / "s.json")
    assert cmd == ["tool", "sub", "cmd", "--mask", str((tmp_path / "m.png").resolve()),
                   "--surface", str((tmp_path / "s.json").resolve())]


def test_pipeline_command_resolves_relative_paths_against_the_working_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cmd = pipeline_command({"cwd": "unused", "argv": ["tool"]}, "m.png", Path("sub") / "s.json")
    assert cmd[2] == str((tmp_path / "m.png").resolve())
    assert cmd[4] == str((tmp_path / "sub" / "s.json").resolve())
    assert Path(cmd[2]).is_absolute() and Path(cmd[4]).is_absolute()


def test_run_pipeline_runs_in_the_configured_directory_with_absolute_inputs(tmp_path, monkeypatch, fake):
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    out = json.loads(run_pipeline(fake("report"), "m.png", "s.json"))
    assert Path(out["cwd"]).resolve() == fake.dir.resolve()
    assert out["args"] == ["--mask", str((work / "m.png").resolve()), "--surface", str((work / "s.json").resolve())]


def test_run_pipeline_without_the_executable_raises_runtime_error(tmp_path):
    cfg = {"cwd": str(tmp_path), "argv": ["no-such-executable-for-crackfem-tests", "run"]}
    with pytest.raises(RuntimeError, match="not found on PATH"):
        run_pipeline(cfg, tmp_path / "m.png", tmp_path / "s.json")


def test_run_pipeline_with_a_missing_directory_raises_file_not_found(tmp_path):
    cfg = {"cwd": str(tmp_path / "no_such_dir"), "argv": [sys.executable, "-V"]}
    with pytest.raises(FileNotFoundError, match="pipeline directory"):
        run_pipeline(cfg, tmp_path / "m.png", tmp_path / "s.json")


def test_failing_pipeline_raises_runtime_error_with_stderr(tmp_path, fake):
    with pytest.raises(RuntimeError, match=r"exit code 3.*boom: cannot read the surface"):
        run_pipeline(fake("fail"), tmp_path / "m.png", tmp_path / "s.json")


def test_pipeline_timeout_raises_runtime_error(tmp_path, fake):
    with pytest.raises(RuntimeError, match="timed out after 1 s"):
        run_pipeline(fake("sleep"), tmp_path / "m.png", tmp_path / "s.json", timeout=1)


@needs_pipeline
def test_live_pipeline_failure_raises_runtime_error_with_stderr(tmp_path):
    write_mask(tmp_path / "m.png", width_px=6)
    with pytest.raises(RuntimeError, match=r"exit code 1.*ENOENT"):
        run_pipeline(load_pipeline_config(), tmp_path / "m.png", tmp_path / "no_such_surface.json")


@needs_pipeline
def test_live_pipeline_accepts_mask_and_surface_paths_relative_to_the_cwd(tmp_path, monkeypatch):
    config = load_pipeline_config()
    write_mask(tmp_path / "m.png", width_px=6)
    write_surface(tmp_path / "s.json", 5e-5)
    monkeypatch.chdir(tmp_path)
    out = parse_pipeline_output(run_pipeline(config, "m.png", "s.json"))
    assert governing_width(out)["max_width_mm"] == pytest.approx(0.30)


# ---- no local path in a saved report ---------------------------------------------------------

HYGIENE_CONFIG = {"cwd": "D:\\somewhere\\tool_dir",
                  "argv": ["tool", "D:\\somewhere\\bin\\entry.js", "pkg/entry.js", "measure"]}


@pytest.mark.parametrize("leak", [
    "D:\\somewhere\\tool_dir",                  # the configured directory, as written
    "D:/somewhere/tool_dir/x",                  # ... with forward slashes
    "D:\\\\somewhere\\\\tool_dir",              # ... JSON-escaped
    "D:\\somewhere\\bin\\entry.js",             # a path-like argv entry
    "pkg/entry.js",                             # a relative path-like argv entry
    "E:\\unrelated\\place",                     # any drive path
    str(ROOT), ROOT.as_posix(), json.dumps(str(ROOT))[1:-1],
], ids=["cwd", "cwd-forward-slashes", "cwd-json-escaped", "argv-absolute", "argv-relative", "other-drive",
        "project", "project-forward-slashes", "project-json-escaped"])
def test_hygiene_check_refuses_a_report_with_a_local_path(leak):
    with pytest.raises(RuntimeError, match="not writing it"):
        run_claimbridge.assert_no_local_paths('{"raw": "GSD: 0.05 mm/px ' + leak + '"}', HYGIENE_CONFIG)


def test_hygiene_check_accepts_a_clean_report_and_ignores_plain_argv_words():
    # "measure" is an argv word, not a path: it must not block "min measurable width"
    text = (FIX / "spike_ok.txt").read_text(encoding="utf-8")
    run_claimbridge.assert_no_local_paths(json.dumps({"raw": text}, ensure_ascii=False), HYGIENE_CONFIG)
