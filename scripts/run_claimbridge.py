"""Run the thesis pipeline on every synthetic case and save the parsed report.

The pipeline is run as configured in pipeline.local.json (or the file named by
CRACK_PIPELINE_CONFIG); see pipeline.example.json. Only the report lines the parser
understands are kept in the saved raw text (the Mask:/Surface: lines carry machine-specific
paths), and the script refuses to write a record that still contains the configured
directory or any configured path, this project's path or any Windows drive path.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem.claimbridge import (load_pipeline_config, parse_pipeline_output,  # noqa: E402
                                  recognised_lines, run_pipeline)

DRIVE_PATH = re.compile(r"(?<![A-Za-z])[A-Za-z]:[\\/]")


def _forms(s: str) -> set[str]:
    """A path string as written, with either slash style, and JSON-escaped."""
    plain = {s, s.replace("\\", "/"), s.replace("/", "\\")}
    return plain | {json.dumps(f)[1:-1] for f in plain}


def _path_like(s: str) -> bool:
    return "/" in s or "\\" in s or bool(DRIVE_PATH.match(s))


def assert_no_local_paths(text: str, config: dict) -> None:
    """Refuse to save text that carries a machine-specific path (the repo is public): the
    configured directory, any path-like entry of the configured command, this project's path
    or any Windows drive path."""
    secrets = [config["cwd"], str(ROOT)] + [a for a in config["argv"] if _path_like(a)]
    for s in secrets:
        if any(form in text for form in _forms(s)):
            raise RuntimeError("a local path leaked into the saved report; not writing it")
    if DRIVE_PATH.search(text.replace("\\\\", "\\")):
        raise RuntimeError("a Windows drive path leaked into the saved report; not writing it")


def main() -> None:
    try:
        config = load_pipeline_config()
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"pipeline configuration problem: {e}")
    data = ROOT / "data" / "claimbridge"
    cases = json.loads((data / "cases.json").read_text(encoding="utf-8"))["cases"]
    for c in cases:
        text = run_pipeline(config, data / f"{c['name']}_mask.png", data / f"{c['name']}_surface.json")
        raw = "\n".join(recognised_lines(text))
        record = {"case": c, "parsed": parse_pipeline_output(text), "raw": raw}
        out = json.dumps(record, indent=2, ensure_ascii=False)
        assert_no_local_paths(out, config)
        (data / f"{c['name']}_spike.json").write_text(out, encoding="utf-8", newline="\n")
        print(c["name"], "->", [(k["gate"], k["max_width_mm"]) for k in record["parsed"]["cracks"]])


if __name__ == "__main__":
    main()
