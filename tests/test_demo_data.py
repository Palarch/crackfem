"""docs/demo_data.js: the browser demo interpolates between solved points; at the report's cases
it must give the report's numbers. The interpolation here is the page's (linear, w ascending)."""
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "demo_data.js"
PREFIX = "window.DEMO_DATA = "


def load(name: str) -> dict:
    return json.loads((ROOT / "results" / name).read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def demo() -> dict:
    text = DATA.read_text(encoding="utf-8")
    body = text[text.index(PREFIX) + len(PREFIX):].strip().rstrip(";")
    return json.loads(body)["models"]


def eta_at(model: dict, w_mm: float) -> float:
    rows = model["rows"][1:]                           # rows[0] is eta = 1; then w rises as eta falls
    return float(np.interp(w_mm, [r["w_mm"] for r in rows], [r["eta"] for r in rows]))


def column_at(model: dict, eta: float, key: str) -> float:
    rows = model["rows"][1:][::-1]                     # eta ascending for np.interp
    return float(np.interp(eta, [r["eta"] for r in rows], [r[key] for r in rows]))


def test_rows_start_uncracked_and_widen_as_eta_falls(demo):
    for m in demo.values():
        rows = m["rows"]
        assert rows[0]["eta"] == 1.0 and rows[0]["w_mm"] == 0.0
        w = [r["w_mm"] for r in rows[1:]]
        assert all(b > a for a, b in zip(w, w[1:]))
        assert len(rows[0]["ei"]) == len(m["x_mid"]) == m["n_el"]


def test_measured_cases_give_the_reported_eta(demo):
    cases = {"single": load("stage1.json")["cases"]}
    s2 = load("stage2.json")
    cases.update({g: s2[g]["cases"] for g in ("G1", "G3")})
    checked = 0
    for name, cs in cases.items():
        for c in cs:
            if c["status"] == "measured":
                assert eta_at(demo[name], c["w_meas_mm"]) == pytest.approx(c["eta_point"], abs=2e-3), (name, c["name"])
                checked += 1
    assert checked >= 9


def test_refused_cases_give_the_reported_lower_bound(demo):
    for c in load("stage1.json")["cases"]:
        if c["status"] == "refused" and c.get("eta_lower_bound_from_gate") is not None:
            gate_mm = 3.0 * c["gsd_mm"]
            assert eta_at(demo["single"], gate_mm) == pytest.approx(c["eta_lower_bound_from_gate"], abs=2e-3), c["name"]


def test_deflection_matches_the_saved_width_response_rows(demo):
    pc = load("proposal_checks.json")
    saved = {"single": pc["stage1"]["rows"], "G1": pc["grillage"]["G1 exterior"]["rows"],
             "G3": pc["grillage"]["G3 interior"]["rows"]}
    for name, rows in saved.items():
        for r in rows[1:]:
            got = column_at(demo[name], r["eta"], "defl_mm")
            assert got == pytest.approx(r["deflection_truck_mm"], rel=0.02, abs=0.2), (name, r["eta"])
