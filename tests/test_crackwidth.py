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


# Hand numbers for PC: rho = 0.6 * 2.6e-3 / 0.13 = 0.012
RHO = 0.012
SR = 3.4 * 0.045 + 1.6 * 0.5 * 0.425 * 12.7e-3 / RHO              # 0.512833 m
C_HAND = 0.4 * 3.5e6 / RHO * (1 + 195 / 35 * RHO)                  # 124.47 MPa


def test_hand_numbers_for_the_effective_ratio_and_both_branches():
    assert PC.rho(AC) == pytest.approx(RHO, rel=1e-12)
    assert PC.width(450e6, AC) == pytest.approx(SR * (450e6 - C_HAND) / 195e9, rel=1e-12)   # main branch
    assert PC.width(200e6, AC) == pytest.approx(SR * 0.6 * 200e6 / 195e9, rel=1e-12)        # 0.6 floor


@pytest.mark.parametrize("f", [0.1, 1.0, 2.4, 2.6, 4.0])
def test_round_trip_straddles_the_branch_boundary(f):
    ds = f * PC._offset(AC)
    assert PC.stress_from_width(PC.width(ds, AC), AC) == pytest.approx(ds, rel=1e-12)


def test_published_ec2_example_intermediates():
    ex = FIX["ec2_crack_example"]
    if ex is None:
        pytest.skip("no open EC2 worked example located (docs/sources.md S21)")
    i = ex["inputs"]
    c = EC2Crack(i["k1"], i["k2"], i["k3"], i["k4"], i["kt"], i["xi"], i["Es"], i["Ecm"],
                 i["fct_eff"], i["cover"], i["phi"], i["As"])
    # printed sigma_s,cr = kt fct,eff / rho (1 + alpha_e rho) = 57.08 MPa
    assert c._offset(i["Ac_eff"]) / 1e6 == pytest.approx(57.08, abs=0.005)
    # printed h = 600 mm, d = 548 mm, x = 237.8 mm -> hc,ef / h = 0.2012
    hc = effective_tension_height(0.6, 0.6 - 0.548, 0.6 - 0.2378)
    assert hc / 0.6 == pytest.approx(ex["printed_intermediates"]["hc_ef_over_h"], abs=5e-5)


def test_width_returns_zero_for_no_tension_before_touching_the_area():
    assert PC.width(0.0, 0.0) == 0.0
    assert PC.width(-50e6, -1.0) == 0.0


@pytest.mark.parametrize("Ac_eff", [0.0, -0.1])
def test_rho_rejects_a_non_positive_effective_area(Ac_eff):
    with pytest.raises(ValueError):
        PC.rho(Ac_eff)


@pytest.mark.parametrize("y_p, y_na", [(0.1, 0.0), (0.1, -0.2), (0.0, 0.9), (-0.05, 0.9)])
def test_effective_tension_height_rejects_non_positive_depths(y_p, y_na):
    with pytest.raises(ValueError):
        effective_tension_height(h=1.5, y_p=y_p, y_na=y_na)
