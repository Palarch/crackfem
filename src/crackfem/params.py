"""Every numeric input of the model. SI base units (N, m, Pa) unless the field name
ends in `_in` or `_ft`: those hold a source formula's native US units and are converted
once where they are used (geometry.py uses INCH, loads.py uses FT).

Each field of Inputs has a source id in SOURCES; every id is a row of
docs/sources.md. Rows tagged [stand-in] are not the real CPK72 girder.
Every value was checked against an opened source in the B1 source pack (2026-10-07).
"""
from dataclasses import dataclass, fields

KIP = 4448.2216152605       # N per kip
FT = 0.3048                 # m per ft
INCH = 0.0254               # m per inch
KSI = 6.894757293168361e6   # Pa per ksi


@dataclass(frozen=True)
class Inputs:
    # span and deck
    span: float = 24.50                 # m
    girder_spacing: float = 2.20        # m
    slab_thickness: float = 0.20        # m
    wearing_thickness: float = 0.05     # m
    n_girders: int = 5                  # count
    n_barriers: int = 2                 # count
    barrier_weight: float = 0.326 * KIP / FT  # N/m per barrier (326 plf, about 4.76 kN/m)
    # girder outline, bottom to top: (height_in, width_bottom_in, width_top_in), inches
    girder_outline_in: tuple[tuple[float, float, float], ...] = (
        (8.0, 26.0, 26.0), (9.0, 26.0, 8.0), (23.0, 8.0, 8.0),
        (6.0, 8.0, 20.0), (8.0, 20.0, 20.0))
    # strands
    strand_diameter: float = 12.7e-3    # m
    strand_area: float = 98.7e-6        # m^2 per strand
    fpu: float = 1860e6                 # Pa
    Ep: float = 195e9                   # Pa
    fpj_ratio: float = 0.75             # jacking stress / fpu, dimensionless
    lump_sum_loss: float = 45.0 * KSI   # Pa
    strand_centroid: float = 0.10       # m above the girder soffit
    strand_cover: float = 0.045         # m, soffit to the bottom strand's surface (EC2 c)
    # concrete
    fck_girder: float = 40e6            # Pa
    fck_slab: float = 28e6              # Pa
    unit_weight_concrete: float = 0.150 * KIP / FT**3  # N/m^3 (150 pcf, about 23.56 kN/m^3)
    unit_weight_wearing: float = 0.140 * KIP / FT**3   # N/m^3 (140 pcf, about 21.99 kN/m^3)
    # EC2 crack width (7.3.4) and deflection (7.4.3)
    k1: float = 1.6                     # EC2, dimensionless
    k2: float = 0.5                     # EC2, dimensionless
    k3: float = 3.4                     # EC2, dimensionless
    k4: float = 0.425                   # EC2, dimensionless
    kt: float = 0.4                     # EC2, dimensionless
    bond_xi: float = 0.6                # EC2, dimensionless
    beta: float = 0.5                   # EC2, dimensionless
    # live load
    axle_loads: tuple[float, ...] = (8.0 * KIP, 32.0 * KIP, 32.0 * KIP)  # N
    axle_spacings: tuple[float, ...] = (14.0 * FT, 14.0 * FT)            # m
    impact_numerator_ft: float = 50.0     # ft
    impact_denominator_ft: float = 125.0  # ft
    impact_cap: float = 0.30              # dimensionless
    df_divisor_ft: float = 5.5            # ft
    live_load_factor: float = 1.30            # DOH design truck = HS20-44 x 1.30 (S12)
    deflection_limit_ratio: float = 800.0  # span / limit
    # grillage (Stage 2)
    overhang_width_ft: float = 3.0         # ft, exterior girder centreline to deck edge
    barrier_width_ft: float = 1.0          # ft, deck edge to the curb face of the barrier
    wheel_spacing_ft: float = 6.0          # ft, transverse distance between the two wheel lines
    truck_clearance_ft: float = 10.0       # ft, clearance width of one truck (wheels 2 ft inside)
    wheel_curb_clearance_ft: float = 2.0   # ft, smallest distance from a wheel line to the curb face
    lane_width_ft: float = 12.0            # ft, design lane width (number of lanes = floor(w / 12))
    poisson_concrete: float = 0.2          # uncracked concrete, dimensionless
    diaphragm_depth: float = 1.10          # m, diaphragm web depth below the slab
    diaphragm_thickness: float = 0.25      # m, diaphragm web thickness


SOURCES = {
    "span": "S01",
    "girder_spacing": "S02", "slab_thickness": "S02", "wearing_thickness": "S02",
    "n_girders": "S02", "n_barriers": "S02", "barrier_weight": "S02",
    "girder_outline_in": "S03",
    "strand_diameter": "S04", "strand_area": "S04", "fpu": "S04",
    "Ep": "S05",
    "fpj_ratio": "S06", "lump_sum_loss": "S06",
    "strand_centroid": "S15", "strand_cover": "S15",
    "fck_girder": "S19", "fck_slab": "S19",
    "unit_weight_concrete": "S11", "unit_weight_wearing": "S11",
    "k1": "S08", "k2": "S08", "k3": "S08", "k4": "S08", "kt": "S08",
    "bond_xi": "S09",
    "beta": "S10",
    "axle_loads": "S12", "axle_spacings": "S12",
    "impact_numerator_ft": "S13", "impact_denominator_ft": "S13", "impact_cap": "S13",
    "df_divisor_ft": "S14",
    "live_load_factor": "S12",
    "deflection_limit_ratio": "S18",
    "overhang_width_ft": "S25", "barrier_width_ft": "S25",
    "wheel_spacing_ft": "S26", "truck_clearance_ft": "S26", "wheel_curb_clearance_ft": "S26",
    "lane_width_ft": "S26",
    "poisson_concrete": "S24",
    "diaphragm_depth": "S27", "diaphragm_thickness": "S27",
}


def field_names() -> set[str]:
    return {f.name for f in fields(Inputs)}
