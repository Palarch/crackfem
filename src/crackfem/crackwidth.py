"""EC2 (EN 1992-1-1) section 7.3.4 crack width for pretensioned members
(EC2 allows delta sigma_p to replace sigma_s) (S08, S09).

  wk = sr,max * (eps_sm - eps_cm)                                           (7.8)
  eps_sm - eps_cm = [ds - kt fct,eff / rho (1 + alpha_e rho)] / Ep >= 0.6 ds / Ep  (7.9)
  rho = rho_p,eff = xi1^2 A_p / Ac,eff, with xi1^2 = xi for prestressing steel only (7.10, 7.5)
  sr,max = k3 c + k1 k2 k4 phi / rho                                         (7.11)
ds = strand stress change from the state of zero concrete strain at the strand level
(delta sigma_p); it replaces sigma_s for prestressed members. For ordinary bars set
xi = 1 and pass the bar area as A_p (used by the published-example test).
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class EC2Crack:
    k1: float
    k2: float
    k3: float
    k4: float
    kt: float
    xi: float
    Ep: float
    Ecm: float
    fct_eff: float
    cover: float
    phi: float
    A_p: float

    def rho(self, Ac_eff: float) -> float:
        if Ac_eff <= 0.0:
            raise ValueError(f"Ac_eff must be positive, got {Ac_eff}")
        return self.xi * self.A_p / Ac_eff

    def sr_max(self, Ac_eff: float) -> float:
        return self.k3 * self.cover + self.k1 * self.k2 * self.k4 * self.phi / self.rho(Ac_eff)

    def _offset(self, Ac_eff: float) -> float:
        r = self.rho(Ac_eff)
        return self.kt * self.fct_eff / r * (1.0 + self.Ep / self.Ecm * r)

    def strain_difference(self, ds: float, Ac_eff: float) -> float:
        if ds <= 0.0:
            return 0.0
        return max((ds - self._offset(Ac_eff)) / self.Ep, 0.6 * ds / self.Ep)

    def width(self, ds: float, Ac_eff: float) -> float:
        if ds <= 0.0:
            return 0.0
        return self.sr_max(Ac_eff) * self.strain_difference(ds, Ac_eff)

    def stress_from_width(self, w: float, Ac_eff: float) -> float:
        """Closed-form inverse of width(); the 0.6 floor governs while ds <= 2.5 * offset."""
        if w <= 0.0:
            return 0.0
        sr = self.sr_max(Ac_eff)
        C = self._offset(Ac_eff)
        if w <= sr * 1.5 * C / self.Ep:
            return w * self.Ep / (0.6 * sr)
        return w * self.Ep / sr + C


def effective_tension_height(h: float, y_p: float, y_na: float) -> float:
    """hc,ef = min(2.5 (h - d), (h - x) / 3, h / 2), with h - d = y_p and h - x = y_na."""
    if y_p <= 0.0 or y_na <= 0.0:
        raise ValueError(f"y_p and y_na must be positive, got y_p={y_p}, y_na={y_na}")
    return min(2.5 * y_p, y_na / 3.0, h / 2.0)
