"""Dead load per girder and the HS20-44 design truck (S02, S11, S12, S13, S14).
Magnitudes are positive here; the caller applies them downward."""
import numpy as np

from .params import FT, Inputs


def stage_a_line(inp: Inputs, girder_area: float) -> float:
    """Girder + wet slab per girder, carried by the bare girder (stage A)."""
    return inp.unit_weight_concrete * (girder_area + inp.slab_thickness * inp.girder_spacing)


def stage_b_line(inp: Inputs) -> float:
    """Wearing surface + barrier share per girder, carried by the composite girder (stage B)."""
    return (inp.unit_weight_wearing * inp.wearing_thickness * inp.girder_spacing
            + inp.n_barriers * inp.barrier_weight / inp.n_girders)


def dead_load_line(inp: Inputs, girder_area: float) -> float:
    return stage_a_line(inp, girder_area) + stage_b_line(inp)


def impact_factor(inp: Inputs) -> float:
    return min(inp.impact_numerator_ft / (inp.span / FT + inp.impact_denominator_ft), inp.impact_cap)


def truck_fraction(inp: Inputs) -> float:
    """Fraction of one truck per girder: (S / 5.5) wheel lines divided by 2 wheel lines."""
    return (inp.girder_spacing / FT) / inp.df_divisor_ft / 2.0


def axle_offsets(inp: Inputs, reverse: bool = False) -> list[tuple[float, float]]:
    offs = np.concatenate([[0.0], np.cumsum(inp.axle_spacings)])
    loads = list(inp.axle_loads)
    if reverse:
        offs, loads = offs[-1] - offs[::-1], loads[::-1]
    return list(zip(offs.tolist(), loads))


def truck_position_max_midspan(inp: Inputs, step: float = 1e-3):
    """Axle positions on the span that maximise the midspan moment, and that moment."""
    L = inp.span
    best_pos, best_m = None, -1.0
    for rev in (False, True):
        axles = axle_offsets(inp, rev)
        offs = np.array([o for o, _ in axles])
        loads = np.array([P for _, P in axles])
        # M(x0) is piecewise linear; it can only peak where an axle sits at midspan or a support,
        # which the uniform grid misses by up to step/2, so add those placements exactly.
        kinks = np.concatenate([L / 2 - offs, -offs, L - offs])
        # the maximum lies at a kink (an axle at midspan or at a support); the 1 mm grid is a belt-and-braces check
        x0 = np.unique(np.concatenate([np.arange(-offs[-1], L + step / 2, step), kinks]))
        X = x0[:, None] + offs[None, :]
        on = (X >= 0.0) & (X <= L)
        il = np.where(X <= L / 2, X / 2, (L - X) / 2)
        m = (np.where(on, il, 0.0) * loads[None, :]).sum(axis=1)
        i = int(np.argmax(m))
        if m[i] > best_m:
            best_m = float(m[i])
            best_pos = [(float(X[i, k]), float(loads[k])) for k in range(loads.size) if on[i, k]]
    return best_pos, best_m


def max_truck_moment_at(inp: Inputs, x: float, step: float = 1e-3) -> float:
    """Unscaled maximum sagging moment (N·m) at section x of the simply supported span over all
    truck positions and both directions (one full truck: not scaled by truck fraction, impact or
    live_load_factor). Influence line IL(s) = s (L - x) / L for s <= x, x (L - s) / L for s > x."""
    L = inp.span
    if not 0.0 <= x <= L:
        raise ValueError(f"x = {x} is outside the span")
    best = 0.0
    for rev in (False, True):
        axles = axle_offsets(inp, rev)
        offs = np.array([o for o, _ in axles])
        loads = np.array([P for _, P in axles])
        # piecewise linear in the truck position: add the exact kinks (an axle at x or at a support)
        kinks = np.concatenate([x - offs, -offs, L - offs])
        x0 = np.unique(np.concatenate([np.arange(-offs[-1], L + step / 2, step), kinks]))
        S = x0[:, None] + offs[None, :]
        on = (S >= 0.0) & (S <= L)
        il = np.where(S <= x, S * (L - x) / L, x * (L - S) / L)
        m = (np.where(on, il, 0.0) * loads[None, :]).sum(axis=1)
        best = max(best, float(m.max()))
    return best


def service_truck_loads(inp: Inputs) -> list[tuple[float, float]]:
    """Axle positions and loads for the service case: already scaled by truck fraction x (1 + impact)
    x live_load_factor (do not scale again). `truck_position_max_midspan` returns the unscaled moment."""
    pos, _ = truck_position_max_midspan(inp)
    s = truck_fraction(inp) * (1.0 + impact_factor(inp)) * inp.live_load_factor
    return [(x, s * P) for x, P in pos]


def midspan_moment_statics(L: float, w: float, points) -> float:
    il = lambda x: x / 2 if x <= L / 2 else (L - x) / 2
    return w * L * L / 8.0 + sum(P * il(x) for x, P in points)
