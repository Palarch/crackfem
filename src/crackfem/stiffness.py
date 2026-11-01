"""Secant flexural stiffness of one element from its moment (S10 main, S20 cross-check).

The FE model carries only the loads applied after a reference state with total moment M_ref
(0 = prestressed and unloaded; the stage-A moment for a staged girder), so
EI_sec = (M - M_ref) / (kappa(M) - kappa(M_ref)), with M and M_ref total moments.
No load increment (M <= M_ref, i.e. unloading or hogging) keeps the uncracked stiffness.
"""
from .section import PCSection, solve_state, uncracked_EI


def _kappa_uncracked(sec: PCSection, eps_p0: float, M: float) -> float:
    return solve_state(sec, eps_p0, M, cracked=False).kappa


def ec2_kappa(sec: PCSection, eps_p0: float, M: float, Mcr: float, beta: float) -> float:
    """EC2 7.4.3: kappa = zeta kappa_II + (1 - zeta) kappa_I, zeta = 1 - beta (Mcr / M)^2."""
    Mcr = max(Mcr, 0.0)
    kI = _kappa_uncracked(sec, eps_p0, M)
    if M <= 0.0 or M <= Mcr:
        return kI
    kII = solve_state(sec, eps_p0, M, cracked=True).kappa
    zeta = 1.0 - beta * (Mcr / M) ** 2
    return zeta * kII + (1.0 - zeta) * kI


def branson_kappa(sec: PCSection, eps_p0: float, M: float, Mcr: float) -> float:
    """Branson-type interpolation on the secant stiffness from kappa(0), used only as a
    cross-check (not ACI's own procedure for prestressed members; ACI 318-19 replaced Branson):
    kappa = kappa(0) + M / (r EI_I + (1 - r) EI_cr), r = (Mcr / M)^3."""
    Mcr = max(Mcr, 0.0)
    if M <= 0.0 or M <= Mcr:
        return _kappa_uncracked(sec, eps_p0, M)
    k0 = _kappa_uncracked(sec, eps_p0, 0.0)
    EI_cr = M / (solve_state(sec, eps_p0, M, cracked=True).kappa - k0)
    r = (Mcr / M) ** 3
    return k0 + M / (r * uncracked_EI(sec) + (1.0 - r) * EI_cr)


def ec2_secant_EI(sec: PCSection, eps_p0: float, M: float, Mcr: float, beta: float,
                  M_ref: float = 0.0) -> float:
    if M <= M_ref or M <= max(Mcr, 0.0):
        return uncracked_EI(sec)
    return (M - M_ref) / (ec2_kappa(sec, eps_p0, M, Mcr, beta) - ec2_kappa(sec, eps_p0, M_ref, Mcr, beta))


def branson_secant_EI(sec: PCSection, eps_p0: float, M: float, Mcr: float, M_ref: float = 0.0) -> float:
    if M <= M_ref or M <= max(Mcr, 0.0):
        return uncracked_EI(sec)
    return (M - M_ref) / (branson_kappa(sec, eps_p0, M, Mcr) - branson_kappa(sec, eps_p0, M_ref, Mcr))
