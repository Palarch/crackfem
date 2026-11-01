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


@pytest.mark.parametrize(
    "y_cut_in, expected_in2",
    [
        (-1, 0.0),
        (0, 0.0),
        (8, 26 * 8),                                # bottom flange only
        (17, 26 * 8 + 9 * (26 + 8) / 2),            # + whole bottom haunch = 361
        (30, 361 + 8 * 13),                         # 13 in into the web = 465
        (40, 361 + 8 * 23),                         # top of web = 545
        (46, 545 + 6 * (8 + 20) / 2),               # + whole top haunch = 629
        (54, 629 + 20 * 8),                         # whole girder = 789
        (60, 789),                                  # cut above the girder, girder_only default
    ],
)
def test_area_below_cut_boundaries(y_cut_in, expected_in2):
    assert area_below(type_iv(), y_cut_in * INCH) / INCH**2 == pytest.approx(expected_in2)


def test_area_below_girder_only_flag_with_slab():
    g = type_iv()
    A_girder = gross_properties(g)[0]
    parts = g + [slab(top_of(g), 0.2, 2.2, E=0.9, fct=1.0)]
    top = top_of(parts)
    assert area_below(parts, top, girder_only=False) == pytest.approx(A_girder + 0.2 * 2.2)
    assert area_below(parts, top) == pytest.approx(A_girder)
