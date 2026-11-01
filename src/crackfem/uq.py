"""Monte Carlo on the measured crack width, and one-at-a-time sensitivity (tornado).

ClaimBridge reports a width w with a band +/- u, where u = 2 GSD (a pixel-resolution band).
We read it as a coverage factor k = 2, i.e. sigma = u / 2 (+/- u covers 95.45 % under a
Gaussian assumption; the exact 95 % factor 1.96 would widen the band by 2 %).

The measured width is drawn from the normal distribution truncated at w > 0, by inverse-CDF
sampling (no rejection loop, so the same seed gives common random numbers across cases).
The band reported for eta is the 5th-95th percentile (a 90 % interval) of the eta values
from measurement noise only, conditional on the model and its parameters. Percentiles are
computed on the samples whose inversion succeeded ("ok"); the other statuses are counted
and reported, and the result is flagged as censored when more than 5 % are dropped.
"""
import math
from dataclasses import dataclass, replace
from typing import Callable

import numpy as np
from scipy.special import ndtr, ndtri

from .inverse import estimate_eta
from .params import Inputs

DEFAULT_SEED = 20261101
CENSORED_FRACTION = 0.05


def _check_width_inputs(w_meas: float, u95: float) -> None:
    if not (math.isfinite(w_meas) and math.isfinite(u95)):
        raise ValueError("w_meas and u95 must be finite")
    if w_meas <= 0.0:
        raise ValueError("w_meas must be > 0")
    if u95 <= 0.0:
        raise ValueError("u95 must be > 0 (u95 == 0 means no uncertainty: use monte_carlo_eta)")


def below_zero_fraction(w_meas: float, u95: float) -> float:
    """Probability mass of N(w_meas, (u95/2)^2) below w = 0, i.e. what truncation discards."""
    _check_width_inputs(w_meas, u95)
    return float(ndtr(-w_meas / (u95 / 2.0)))


def sample_widths(w_meas: float, u95: float, n: int, seed: int = DEFAULT_SEED) -> np.ndarray:
    """n draws from N(w_meas, (u95/2)^2) truncated at w > 0 (inverse CDF)."""
    _check_width_inputs(w_meas, u95)
    sigma = u95 / 2.0
    lo = ndtr(-w_meas / sigma)
    u = np.random.default_rng(seed).uniform(lo, 1.0, size=n)
    w = w_meas + sigma * ndtri(u)
    return np.maximum(w, np.nextafter(0.0, 1.0))


def summarize(values) -> dict:
    v = np.asarray(values, dtype=float)
    p05, p50, p95 = np.percentile(v, [5, 50, 95])
    return {"n": int(v.size), "mean": float(v.mean()), "p05": float(p05), "p50": float(p50), "p95": float(p95)}


@dataclass(frozen=True)
class MonteCarlo:
    widths: np.ndarray
    etas: np.ndarray            # only samples with status "ok"
    status_counts: dict
    status: np.ndarray          # status string of every sample, aligned with widths
    below_zero_fraction: float  # normal mass below w = 0 discarded by the truncation

    def summary(self) -> dict | None:
        """Percentiles of the ok samples plus how many were dropped; None if no sample is ok."""
        if not self.etas.size:
            return None
        n_total = int(self.status.size)
        n_dropped = n_total - int(self.etas.size)
        out = summarize(self.etas)
        out.update(n_total=n_total, n_dropped=n_dropped, status_counts=dict(self.status_counts),
                   censored=n_dropped / n_total > CENSORED_FRACTION)
        return out


def monte_carlo_eta(w_pred: Callable[[float], float], w_meas: float, u95: float, n: int = 2000,
                    seed: int = DEFAULT_SEED, **estimate_kw) -> MonteCarlo:
    if u95 == 0.0:
        est = estimate_eta(w_pred, float(w_meas), **estimate_kw)   # one solve fills all n
        widths = np.full(n, float(w_meas))
        status = np.full(n, est.status)
        etas = np.full(n, est.eta) if est.status == "ok" else np.empty(0)
        bz = 0.0
    else:
        widths = sample_widths(w_meas, u95, n, seed)
        bz = below_zero_fraction(w_meas, u95)
        etas_list, status_list = [], []
        for w in widths:
            est = estimate_eta(w_pred, float(w), **estimate_kw)
            status_list.append(est.status)
            if est.status == "ok":
                etas_list.append(est.eta)
        status = np.array(status_list, dtype=str)
        etas = np.array(etas_list)
    names, counts = np.unique(status, return_counts=True)
    status_counts = {str(k): int(c) for k, c in zip(names, counts)}
    return MonteCarlo(widths, etas, status_counts, status, bz)


def tornado(metric: Callable[[Inputs], float], base: Inputs, perturbations: dict) -> list[dict]:
    """perturbations: {field name: (low value, high value)}; rows sorted by swing.

    metric must return a float. swing = max - min of the metric over (low, base, high), so a
    parameter whose effect is non-monotone (extremum at the base value) is not under-ranked.

    If any of the three metric values is not finite (e.g. the inverse failed at one end of the
    range), the swing is NaN and the row carries "incomplete": True, because max/min with a NaN
    would silently return a number that depends on the argument order. Incomplete rows are
    listed last, after the rows sorted by descending finite swing."""
    m0 = metric(base)
    rows = []
    for name, (lo, hi) in perturbations.items():
        m_lo = metric(replace(base, **{name: lo}))
        m_hi = metric(replace(base, **{name: hi}))
        row = {"name": name, "low": lo, "high": hi, "base": m0, "metric_low": m_lo, "metric_high": m_hi}
        if all(math.isfinite(m) for m in (m_lo, m0, m_hi)):
            row["swing"] = max(m_lo, m0, m_hi) - min(m_lo, m0, m_hi)
        else:
            row["swing"] = math.nan
            row["incomplete"] = True
        rows.append(row)
    return sorted(rows, key=lambda r: (r.get("incomplete", False), 0.0 if r.get("incomplete") else -r["swing"]))
