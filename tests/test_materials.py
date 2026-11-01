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
