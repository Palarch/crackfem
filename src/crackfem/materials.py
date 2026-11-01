"""Concrete material formulas, Pa in and out.

EC2 (EN 1992-1-1 Table 3.1, source S07) is the main family. The ACI 318-19 formulas
(source S16: aci_ec, aci_fr) are kept for reference and tests; the Branson cross-check
does not use them (it applies the Branson interpolation to the same EC2 section and Mcr).
"""
import math


def ec2_fcm(fck: float) -> float:
    """Mean compressive strength: fcm = fck + 8 MPa."""
    return fck + 8e6


def ec2_fctm(fck: float) -> float:
    """Mean axial tensile strength: fctm = 0.30 fck^(2/3) (MPa), valid up to C50/60."""
    if fck > 50e6:
        raise ValueError("fctm = 0.30 fck^(2/3) is valid only up to C50/60")
    return 0.30e6 * (fck / 1e6) ** (2.0 / 3.0)


def ec2_ecm(fck: float) -> float:
    """Secant modulus: Ecm = 22 (fcm / 10)^0.3 GPa, fcm in MPa."""
    return 22e9 * (ec2_fcm(fck) / 1e6 / 10.0) ** 0.3


def aci_ec(fc: float) -> float:
    """ACI 318-19 normal-weight concrete: Ec = 4700 sqrt(f'c) MPa."""
    return 4700.0 * math.sqrt(fc / 1e6) * 1e6


def aci_fr(fc: float) -> float:
    """ACI 318-19 modulus of rupture: fr = 0.62 sqrt(f'c) MPa."""
    return 0.62 * math.sqrt(fc / 1e6) * 1e6
