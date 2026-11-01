from crackfem import figures as F


def test_width_response_figure_writes_a_png(tmp_path):
    p = tmp_path / "w.png"
    F.fig_width_response({"girder": ([0.0, 0.1, 0.4], [9.5, 60.0, 80.0]), "G1": ([0.0, 0.2], [8.0, 20.0])},
                         p, ylabel="truck deflection (mm)", hlines={"L/800 = 30.6 mm": 30.6},
                         vlines={"3 px at 0.02 mm/px": 0.06})
    assert p.stat().st_size > 1000


def test_damaged_length_figure_writes_a_png(tmp_path):
    p = tmp_path / "d.png"
    F.fig_damaged_length([10, 20, 40], {"one element": [30.0, 20.0, 15.0], "fixed 4.9 m": [40.0, 40.5, 40.6]},
                         [0.3, 2.45, 24.5], [12.0, 30.0, 72.0], p,
                         hlines={"distributed": 72.0, "intact": 9.5}, marks={"N = 10 element": 2.45})
    assert p.stat().st_size > 1000


def test_width_response_figure_draws_measurement_bands(tmp_path):
    p = tmp_path / "b.png"
    F.fig_width_response({"girder": ([0.0, 0.1, 0.4], [9.5, 60.0, 80.0])}, p, ylabel="truck deflection (mm)",
                         bands=[{"curve": "girder", "w": 0.16, "u": 0.04, "lo": 69.0, "mid": 72.0, "hi": 74.4}])
    assert p.stat().st_size > 1000
