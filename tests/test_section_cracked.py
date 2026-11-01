import json
import math
from pathlib import Path

import numpy as np
import pytest

from crackfem.geometry import Trapezoid, discretize, slab, top_of
from crackfem.section import (PCSection, concrete_stress, delta_sigma_p, resultants, soffit_stress,
                              solve_state, strand_stress)

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


def test_between_top_and_soffit_decompression_the_cracked_state_equals_the_uncracked_state():
    # This section's strand lies below the kern, so prestress alone puts the top in tension:
    # the modes agree only between top decompression (~185 kN·m) and soffit decompression (~503 kN·m).
    sec = pc_section()
    a = solve_state(sec, 4e-3, 300e3, cracked=False)
    b = solve_state(sec, 4e-3, 300e3, cracked=True)
    assert (b.eps0, b.kappa) == pytest.approx((a.eps0, a.kappa), rel=1e-12)


def test_prestress_alone_cracks_the_top_when_the_strand_is_below_the_kern():
    sec = pc_section()
    a = solve_state(sec, 4e-3, 50e3, cracked=False)
    b = solve_state(sec, 4e-3, 50e3, cracked=True)
    assert b.kappa < a.kappa                       # more hogging once the top tension is lost
    top = sec.fib.y > 0.79
    assert concrete_stress(sec, b)[top].max() == 0.0


@pytest.mark.parametrize("cracked", [False, True])
def test_unloaded_section_without_prestress_returns_the_zero_state(cracked):
    st = solve_state(pc_section(), 0.0, 0.0, cracked=cracked)
    assert (st.eps0, st.kappa, st.M) == (0.0, 0.0, 0.0)
    assert st.cracked is cracked


def test_soffit_stress_rejects_a_cracked_state():
    sec = pc_section()
    with pytest.raises(ValueError):
        soffit_stress(sec, solve_state(sec, 4e-3, 700e3, cracked=True))


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
        assert h_tot - st.neutral_axis == pytest.approx(ex["neutral_axis_depth_m"], rel=1e-3)
    if ex["strand_stress_MPa"] is not None:
        assert strand_stress(sec, e0, st) / 1e6 == pytest.approx(ex["strand_stress_MPa"], rel=0.03)
    if ex.get("delta_sigma_p_MPa") is not None:
        assert delta_sigma_p(sec, st) / 1e6 == pytest.approx(ex["delta_sigma_p_MPa"], rel=5e-3)
