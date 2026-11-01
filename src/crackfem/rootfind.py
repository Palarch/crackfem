"""Illinois method: a secant step that always keeps the root bracketed (the 'safeguarded
secant' of the spec). When the same end survives twice, its f value is halved so that
both ends keep moving.

xtol must be well above one ulp of x: the stopping test is |b - a| <= xtol, and a bracket
cannot be narrower than the spacing of the floating-point numbers around the root.
Non-finite f values (NaN, inf) are refused with a ValueError instead of being followed.
"""
import math
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class RootResult:
    x: float
    iterations: int
    history: tuple      # (x, f(x)) for both starting ends and every iterate


def illinois(f: Callable[[float], float], a: float, b: float,
             xtol: float = 1e-10, max_iter: int = 300, *,
             fa: float | None = None, fb: float | None = None) -> RootResult:
    """Root of f in [a, b]. If the caller already knows f(a) / f(b) it can pass them as
    fa / fb to save the two end evaluations; the history still records both ends."""
    if fa is None:
        fa = f(a)
    if fb is None:
        fb = f(b)
    history = [(a, fa), (b, fb)]
    if not (math.isfinite(fa) and math.isfinite(fb)):
        raise ValueError("f(a) or f(b) is not finite")
    if fa == 0.0:
        return RootResult(a, 0, tuple(history))
    if fb == 0.0:
        return RootResult(b, 0, tuple(history))
    if not (fa * fb < 0.0):
        raise ValueError("root is not bracketed: f(a) and f(b) have the same sign")
    side = 0
    for it in range(1, max_iter + 1):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        history.append((c, fc))
        if not math.isfinite(fc):
            raise ValueError(f"f({c!r}) is not finite")
        if fc == 0.0:
            return RootResult(c, it, tuple(history))
        if fc * fb > 0.0:          # root in [a, c]: c replaces b, a survives
            b, fb = c, fc
            if side == -1:
                fa *= 0.5
            side = -1
        else:                      # root in [c, b]: c replaces a, b survives
            a, fa = c, fc
            if side == +1:
                fb *= 0.5
            side = +1
        if abs(b - a) <= xtol:
            return RootResult(0.5 * (a + b), it, tuple(history))
    raise RuntimeError(f"Illinois did not converge in {max_iter} iterations")
