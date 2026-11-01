"""Synthetic scenarios: eta_true -> predicted crack width -> masks at several GSDs.

GSDs: 0.02 and 0.05 mm/px (close-range photo), 0.10 mm/px, 0.38 mm/px (lower end of the
thesis pipeline's field GSD). Widths are rounded to whole pixels, as a camera would.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crackfem.claimbridge import write_mask, write_surface  # noqa: E402
from crackfem.girder1d import Girder1D  # noqa: E402
from crackfem.inverse import eta_decompression  # noqa: E402

GSDS_MM = (0.02, 0.05, 0.10, 0.38)
DROPS = (0.10, 0.20, 0.35)    # eta_true = eta_decompression - drop; 0.35 is measurable at GSD 0.10


def main() -> None:
    g = Girder1D()
    eta_dec = eta_decompression(g.w_pred)
    out = ROOT / "data" / "claimbridge"
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    for k, drop in enumerate(DROPS, start=1):
        eta_true = round(eta_dec - drop, 2)
        w_true_mm = g.w_pred(eta_true) * 1e3
        for gsd in GSDS_MM:
            width_px = max(1, round(w_true_mm / gsd))
            name = f"case{k}_gsd{gsd:.2f}".replace(".", "p")
            write_mask(out / f"{name}_mask.png", width_px)
            write_surface(out / f"{name}_surface.json", gsd * 1e-3)
            cases.append({"name": name, "eta_true": eta_true, "w_true_mm": w_true_mm,
                          "gsd_mm": gsd, "width_px": width_px})
    (out / "cases.json").write_text(json.dumps({"eta_decompression": eta_dec, "cases": cases}, indent=2),
                                    encoding="utf-8", newline="\n")
    for c in cases:
        print(f"{c['name']}: eta_true={c['eta_true']:.2f} w_true={c['w_true_mm']:.3f} mm -> {c['width_px']} px")


if __name__ == "__main__":
    main()
