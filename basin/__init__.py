"""basin: a law finder for tables of measurements.

BACON proposes terms (including fractional powers), sparse regression builds candidate
laws y = N(x)/D(x), and cross-validation decides between them.

    >>> import numpy as np, basin
    >>> law = basin.discover(X, y, ["M1", "T1", "M2", "T2"])

The method itself lives in `basin.core`, unchanged from the benchmarked version; its
module-level switches (`POWER_FREE`, `POWER_DENOMINATORS`, ...) can be changed for one
call with `basin.settings(...)`.
"""
from contextlib import contextmanager

import sympy as sp

from basin import core
from basin.core import PARAMS, propose_powers, term_list
from basin.version import __version__

__all__ = ["discover", "settings", "propose_powers", "term_list", "core", "PARAMS",
           "SWITCHES", "__version__"]

# Settings that `settings()` may change. Their defaults are the benchmark version.
SWITCHES = ("MAX_FREE", "EXHAUSTIVE_TERMS", "SCREEN_KEEP", "POWER_SE", "POWER_FLOOR",
            "POWER_KEEP", "CV_FLOOR", "FOLDS", "REPEATS", "CV_SEED", "POWER_DENOMINATORS",
            "POWER_FREE", "POWER_FREE_DIGITS")


def _symbols(symbols, n):
    if symbols is None:
        return [sp.Symbol(f"x{i}", positive=True) for i in range(n)]
    if len(symbols) != n:
        raise ValueError(f"{len(symbols)} names given for {n} input columns")
    return [s if isinstance(s, sp.Symbol) else sp.Symbol(str(s), positive=True)
            for s in symbols]


def discover(X, y, symbols=None, details=False):
    """Find a law y = N(x)/D(x) for a table of measurements.

    Parameters
    ----------
    X : array-like, shape (n_rows, n_inputs)
        One row per measurement, one column per input.
    y : array-like, shape (n_rows,)
        The measured output.
    symbols : list of str or sympy.Symbol, optional
        Names of the inputs. Strings become positive sympy symbols. Default x0, x1, ...
    details : bool
        If True, also return a dict with the proposed powers, their estimates, the size of
        the term list, the number of candidate laws and the chosen law's CV error.

    Returns
    -------
    sympy.Expr, or (sympy.Expr, dict) if details=True. None (or (None, {})) if no
    candidate law could be fitted.
    """
    import numpy as np

    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    return core.discover(X, y, _symbols(symbols, X.shape[1]), details=details)


@contextmanager
def settings(**switches):
    """Temporarily change the method's switches, e.g.

        with basin.settings(POWER_FREE=True):
            law = basin.discover(X, y, names)

    The previous values are restored on exit. Not thread-safe (the switches are module
    globals, as in the benchmarked code).
    """
    unknown = set(switches) - set(SWITCHES)
    if unknown:
        raise ValueError(f"unknown switch(es): {sorted(unknown)}; known: {list(SWITCHES)}")
    old = {k: getattr(core, k) for k in switches}
    try:
        for k, v in switches.items():
            setattr(core, k, v)
        yield
    finally:
        for k, v in old.items():
            setattr(core, k, v)
