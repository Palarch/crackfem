import math

import pytest
from scipy.optimize import brentq

from crackfem.rootfind import illinois

CASES = [
    (lambda x: x**3 - 2 * x - 5, 2.0, 3.0),
    (lambda x: math.cos(x) - x, 0.0, 1.0),
    (lambda x: 2 * max(0.0, 0.6 - x) - 0.1, 0.0, 1.0),   # flat beyond 0.6, like w_pred(eta)
]


@pytest.mark.parametrize("f,a,b", CASES)
def test_matches_brent(f, a, b):
    assert illinois(f, a, b, xtol=1e-12).x == pytest.approx(brentq(f, a, b, xtol=1e-14), abs=1e-10)


def test_refuses_an_unbracketed_interval():
    with pytest.raises(ValueError):
        illinois(lambda x: x * x + 1, -1.0, 1.0)


def test_records_history_and_counts_iterations():
    r = illinois(lambda x: x - 0.3, 0.0, 1.0)
    assert r.iterations >= 1
    assert len(r.history) == r.iterations + 2


NAN = float("nan")


def test_nan_at_an_end_is_refused():
    with pytest.raises(ValueError):
        illinois(lambda x: NAN if x == 0.0 else x - 0.3, 0.0, 1.0)
    with pytest.raises(ValueError):
        illinois(lambda x: NAN if x == 1.0 else x - 0.3, 0.0, 1.0)


def test_nan_inside_is_refused_at_once():
    calls = []

    def f(x):
        calls.append(x)
        return x - 0.3 if x in (0.0, 1.0) else NAN

    with pytest.raises(ValueError):
        illinois(f, 0.0, 1.0)
    assert len(calls) == 3          # two ends and the first iterate, then it stops


def test_infinite_known_end_values_are_refused():
    with pytest.raises(ValueError):
        illinois(lambda x: x - 0.3, 0.0, 1.0, fa=float("inf"), fb=0.7)


def test_known_end_values_save_the_end_evaluations_but_stay_in_the_history():
    calls = []

    def f(x):
        calls.append(x)
        return x - 0.3

    r = illinois(f, 0.0, 1.0, fa=-0.3, fb=0.7)
    assert r.x == pytest.approx(0.3, abs=1e-10)
    assert 0.0 not in calls and 1.0 not in calls
    assert r.history[:2] == ((0.0, -0.3), (1.0, 0.7))
    assert len(r.history) == r.iterations + 2


def test_a_root_exactly_at_an_end_returns_without_iterating():
    r = illinois(lambda x: x - 1.0, 0.0, 1.0)
    assert r.x == 1.0 and r.iterations == 0
    r = illinois(lambda x: x, 0.0, 1.0)
    assert r.x == 0.0 and r.iterations == 0
