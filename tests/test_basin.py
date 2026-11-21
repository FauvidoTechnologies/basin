from itertools import product

import numpy as np
import pytest
import sympy as sp

import basin


def _black(noise=0.0):
    levels = [[1.0, 2.0, 3.0], [10.0, 50.0, 90.0], [1.0, 2.0, 3.0], [20.0, 60.0, 100.0]]
    X = np.array(list(product(*levels)))
    M1, T1, M2, T2 = X.T
    y = (M1 * T1 + M2 * T2) / (M1 + M2)
    if noise:
        y = y + np.random.default_rng(0).normal(0.0, noise * np.abs(y))
    return X, y


def test_black_law_exact():
    X, y = _black()
    law = basin.discover(X, y, ["M1", "T1", "M2", "T2"])
    M1, T1, M2, T2 = sp.symbols("M1 T1 M2 T2", positive=True)
    f = sp.lambdify((M1, T1, M2, T2), law)
    assert np.allclose(f(*X.T), y, rtol=1e-6)


def test_power_law_proposes_fraction():
    x = np.repeat([1.0, 2.0, 4.0, 8.0, 16.0], 1)
    X = np.array(list(product(x, [1.0, 2.0, 3.0])))
    y = 2.0 * X[:, 0] ** 1.5 / X[:, 1]
    law, info = basin.discover(X, y, details=True)
    assert "3/2" in info["proposed_powers"]["x0"]
    f = sp.lambdify(sp.symbols("x0 x1", positive=True), law)
    assert np.allclose(f(*X.T), y, rtol=1e-6)


def test_hybrid_shim_is_core():
    import hybrid
    assert hybrid is basin.core


def test_settings_restores():
    assert basin.core.POWER_FREE is False
    with basin.settings(POWER_FREE=True, POWER_DENOMINATORS=(1, 2, 3, 4)):
        assert basin.core.POWER_FREE is True
    assert basin.core.POWER_FREE is False
    assert basin.core.POWER_DENOMINATORS == (1, 2, 3)
    with pytest.raises(ValueError):
        with basin.settings(NOT_A_SWITCH=1):
            pass


def test_symbol_count_checked():
    X, y = _black()
    with pytest.raises(ValueError):
        basin.discover(X, y, ["a", "b"])
