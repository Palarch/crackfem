from crackfem import figures as F


def test_every_figure_function_writes_a_png(tmp_path):
    case = {"name": "c", "gsd_mm": 0.05, "eta_true": 0.4, "status": "measured", "w_meas_mm": 0.15,
            "u95_mm": 0.10, "eta_point": 0.41, "eta_mc": {"p05": 0.3, "p50": 0.41, "p95": 0.5}}
    refused = {**case, "name": "r", "gsd_mm": 0.38, "status": "refused", "eta_point": None, "eta_mc": None}
    calls = [
        (F.fig_w_eta, ([0.1, 0.5, 1.0], [0.4, 0.1, 0.0], [case])),
        (F.fig_mkappa, ([0, 1000, 2000], {"uncracked": [0, 1e-4, 2e-4], "cracked": [0, 2e-4, 5e-4]})),
        (F.fig_ei_ratio, ([1.0, 2.0, 3.0], {"eta = 1": [1, 1, 1], "eta = 0.4": [1, 0.6, 1]})),
        (F.fig_deflection, ([0, 1, 2], {"eta = 1": [0, 10, 0], "eta = 0.4": [0, 14, 0]})),
        (F.fig_eta_hist, ([0.3, 0.4, 0.41, 0.5], 0.4)),
        (F.fig_mesh, ([10, 20, 40], {"intact": [5, 5, 5], "damaged": [7, 7.2, 7.25]})),
        (F.fig_cases, ([case, refused],)),
        (F.fig_cases, ([case, {**refused, "name": "u", "status": "unresolved"}],)),
        (F.fig_tornado, ([{"name": "a", "low": 1, "high": 2, "base": 0.4, "metric_low": 0.35, "metric_high": 0.45}],)),
    ]
    for i, (fn, args) in enumerate(calls):
        p = tmp_path / f"f{i}.png"
        fn(*args, p)
        assert p.stat().st_size > 1000, fn.__name__
