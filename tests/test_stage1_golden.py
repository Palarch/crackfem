"""Golden numbers of the Stage-1 fallback (results/stage1.json, 2026-10-12): refactoring for Stage 2
must not change them silently.

Computed fresh from the model (results/ is not read), so the test is fast and fails when the
physics or the design rule moves, not when a result file is stale."""
import pytest

from crackfem.girder1d import Girder1D
from crackfem.inverse import eta_decompression


@pytest.fixture(scope="module")
def g():
    return Girder1D()


def test_design_has_30_strands(g):
    assert g.design.n_strands == 30


def test_decompression_eta(g):
    assert eta_decompression(g.w_pred) == pytest.approx(0.5566, abs=5e-4)


def test_predicted_crack_width_at_eta_036(g):
    assert g.w_pred(0.36) * 1e3 == pytest.approx(0.1612, rel=2e-3)          # mm


def test_intact_truck_deflection(g):
    assert g.indicator(1.0)["deflection_truck_mm"] == pytest.approx(9.506, rel=1e-3)


def test_damaged_truck_deflection_and_stiffness_ratio(g):
    ind = g.indicator(0.36)
    assert ind["deflection_truck_mm"] == pytest.approx(72.12, rel=2e-3)
    assert min(ind["EI_ratio"]) == pytest.approx(0.1217, rel=5e-3)
