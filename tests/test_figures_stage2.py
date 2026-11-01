from crackfem import figures as F


def test_every_stage2_figure_function_writes_a_png(tmp_path):
    band = {"p05": 0.3, "p50": 0.35, "p95": 0.4}
    rows = [{"name": "c1", "eta_gen_1d": 0.36, "bands": {"Stage 1": band, "grillage G1": band, "grillage G3": None}},
            {"name": "c2", "eta_gen_1d": 0.21, "bands": {"Stage 1": band, "grillage G1": band, "grillage G3": band}}]
    calls = [
        lambda p: F.fig_w_eta_models([0.1, 0.5, 1.0], {"a": [0.4, 0.1, 0.0], "b": [0.3, 0.0, 0.0]}, p,
                                     eta_dec={"a": 0.6}),
        lambda p: F.fig_eta_bands(rows, p, title="t"),
        lambda p: F.fig_girder_bars(["G1", "G2"], {"intact": [1.0, 2.0], "damaged": [0.5, 2.6]}, p, "M", title="t"),
        lambda p: F.fig_history({"G1": [0.5, 0.1, 1e-4, 1e-7], "G3": [0.7, 0.0]}, p, title="t"),
    ]
    for i, draw in enumerate(calls):
        p = tmp_path / f"s{i}.png"
        draw(p)
        assert p.stat().st_size > 1000, i


def test_limit_line_and_mesh_label_options_write_a_png(tmp_path):
    p = tmp_path / "bars.png"
    F.fig_girder_bars(["G1", "G2"], {"intact": [9.5, 8.0], "damaged": [21.9, 15.2]}, p, "mm",
                      hline=(30.6, "limit L/800 = 30.6 mm"))
    assert p.stat().st_size > 1000
    q = tmp_path / "mesh.png"
    F.fig_mesh([10, 20, 40], {"a": [1.0, 1.01, 1.02]}, q, ylabel="relative change")
    assert q.stat().st_size > 1000
