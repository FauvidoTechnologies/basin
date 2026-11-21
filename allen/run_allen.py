"""Allen Cell Types: current needed to fire (rheobase) from input resistance, time
constant and the voltage gap to threshold. Hybrid vs PySINDy vs PySR.

Data: Allen Brain Map API, model::EphysFeature (2,336 cells), saved in data/allen/.
    R   = ri                        input resistance, MOhm (or input_resistance_mohm, --rcol)
    tau = tau                       membrane time constant, ms
    gap = threshold_v_long_square - vrest     mV
    I   = threshold_i_long_square   rheobase, pA (long 1 s square pulses)
Known answer: Ohm's law at steady state, I = gap / R, i.e. I[pA] = 1000 * gap[mV] / R[MOhm];
tau drops out (a 1 s step reaches steady state).

One job = one method, one resistance column, one fold:
    python run_allen.py job --method Hybrid --rcol ri --fold all --out results/
fold "all" fits every cell (the reported formula); folds 0-4 fit 80% and are scored on
the other 20% (cells split once, seed 0). Then:
    python run_allen.py score --out results/

Settings are the bacon-vs-pysindy ones, unchanged:
  Hybrid   hybrid.discover, benchmark defaults.
  PySINDy  run_static.sindy_static with FEYNMAN_LIB (integer powers up to +-3 per input,
           total degree <= 5, STLSQ normalized, threshold sweep, sparsest within
           R2_TOL of the best RSS), copied here so this file does not import BACON.
  PySR     2.7.2, run_static.PYSR_SETTINGS + PYSR_FEYNMAN (+ - * / sqrt exp sin cos asin,
           maxsize 40, 100 iterations, deterministic, serial). One search; both the
           "best" (PySR default) and "score" picks are recorded, as in Test 17.
"""
import argparse
import glob
import json
import os
import sys
import time
import warnings
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import sympy as sp  # noqa: E402

NAMES = ["R", "tau", "gap"]
SYMS = {v: sp.Symbol(v, positive=True) for v in NAMES}
N_FOLDS, SPLIT_SEED = 5, 0
PHYSICAL_C = 1000.0          # pA per (mV / MOhm)
COEF_TOL = 0.05              # "exact law": coefficient within 5%, as in the benchmarks


# --- data ------------------------------------------------------------------------------
def load(rcol, data_dir=os.path.join(REPO, "data", "allen")):
    with open(os.path.join(data_dir, "ephys_features.json")) as fh:
        rows = json.load(fh)["msg"]
    rows = sorted(rows, key=lambda r: r["specimen_id"])
    keep = [r for r in rows if all(r.get(k) is not None for k in
            (rcol, "tau", "vrest", "threshold_v_long_square", "threshold_i_long_square"))]
    X = np.array([[r[rcol], r["tau"], r["threshold_v_long_square"] - r["vrest"]]
                  for r in keep], dtype=float)
    y = np.array([r["threshold_i_long_square"] for r in keep], dtype=float)
    ok = np.all(X > 0, axis=1) & (y > 0)
    return X[ok], y[ok]


def folds(n):
    perm = np.random.default_rng(SPLIT_SEED).permutation(n)
    return np.array_split(perm, N_FOLDS)


# --- methods ---------------------------------------------------------------------------
def _syms(names):
    return {v: sp.Symbol(v, positive=True) for v in names}


# Hybrid side runs: named switches in hybrid.py (defaults = the benchmark version).
HYBRID_VARIANTS = {"Hybrid": {},
                   "Hybrid-den6": {"POWER_DENOMINATORS": (1, 2, 3, 4, 5, 6)},
                   "Hybrid-free": {"POWER_FREE": True}}


def run_hybrid(X, y, names=NAMES, variant="Hybrid"):
    import hybrid
    syms = _syms(names)
    old = {k: getattr(hybrid, k) for k in HYBRID_VARIANTS[variant]}
    for k, v in HYBRID_VARIANTS[variant].items():
        setattr(hybrid, k, v)
    try:
        return {variant: hybrid.discover(X, y, [syms[v] for v in names])}
    finally:
        for k, v in old.items():
            setattr(hybrid, k, v)


def _exponents(n, degree, max_abs):
    out = [e for e in product(range(-max_abs, max_abs + 1), repeat=n)
           if sum(abs(v) for v in e) <= degree]
    return sorted(out, key=lambda e: (sum(abs(v) for v in e), [-v for v in e]))


def run_pysindy(X, y, degree=5, max_abs=3, r2_tol=0.01, names=NAMES):
    import pysindy as ps
    from hybrid import THRESHOLDS_NORM          # = sindy_runner.THRESHOLDS_NORM
    exps = _exponents(X.shape[1], degree, max_abs)
    theta = np.column_stack([np.prod(X ** np.array(e, dtype=float), axis=1) for e in exps])
    total_var = np.sum((y - y.mean()) ** 2)
    fits = []
    for thr in THRESHOLDS_NORM:
        opt = ps.STLSQ(threshold=thr, normalize_columns=True, alpha=0.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            opt.fit(theta, y.reshape(-1, 1))
        coef = np.asarray(opt.coef_).ravel()
        k = int(np.count_nonzero(coef))
        if k > 0 and np.all(np.isfinite(coef)):
            fits.append((k, float(np.sum((y - theta @ coef) ** 2)), coef))
    if not fits:
        return {"PySINDy": None}
    best = min(f[1] for f in fits)
    coef = min([f for f in fits if f[1] <= best + r2_tol * total_var],
               key=lambda f: (f[0], f[1]))[2]
    syms = [_syms(names)[v] for v in names]
    terms = [sp.Mul(*[s ** int(p) for s, p in zip(syms, e)]) for e in exps]
    return {"PySINDy": sum((float(c) * t for c, t in zip(coef, terms) if c != 0),
                           sp.Integer(0))}


PYSR_SETTINGS = {"deterministic": True, "parallelism": "serial", "random_state": 0}
PYSR_FEYNMAN = {"binary_operators": ["+", "-", "*", "/"],
                "unary_operators": ["sqrt", "exp", "sin", "cos", "asin"], "maxsize": 40}


def run_pysr(X, y, names=NAMES):
    from pysr import PySRRegressor
    from pysr.export_sympy import assert_valid_sympy_symbol
    from pysr.sr import idx_model_selection
    syms = _syms(names)
    alias = {}
    for v in names:
        try:
            assert_valid_sympy_symbol(v, None)
        except ValueError:
            alias[v] = f"{v}_"
    m = PySRRegressor(**PYSR_SETTINGS, **PYSR_FEYNMAN, verbosity=0, progress=False,
                      temp_equation_file=True, delete_tempfiles=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m.fit(X, y, variable_names=[alias.get(v, v) for v in names])
    eq = m.equations_
    back = {a: v for v, a in alias.items()}
    out, kept = {}, []
    for rule in ("best", "score"):
        e = eq.loc[idx_model_selection(eq, rule), "sympy_format"]
        out[f'PySR "{rule}"'] = e.xreplace({s: syms[back.get(s.name, s.name)]
                                            for s in e.free_symbols})
    for _, row in eq.iterrows():
        kept.append({"complexity": int(row["complexity"]), "loss": float(row["loss"]),
                     "score": float(row["score"]), "equation": str(row["equation"])})
    return out, kept


# --- one job ---------------------------------------------------------------------------
def job(method, rcol, fold, out):
    X, y = load(rcol)
    tr = np.arange(len(y)) if fold == "all" else np.setdiff1d(
        np.arange(len(y)), folds(len(y))[int(fold)])
    t0 = time.time()
    kept = None
    if method.startswith("Hybrid"):
        found = run_hybrid(X[tr], y[tr], variant=method)
    elif method == "PySINDy":
        found = run_pysindy(X[tr], y[tr])
    else:
        found, kept = run_pysr(X[tr], y[tr])
    rec = {"method": method, "rcol": rcol, "fold": fold, "n_train": int(len(tr)),
           "seconds": round(time.time() - t0, 2),
           "formulas": {k: (None if e is None else sp.srepr(e)) for k, e in found.items()},
           "readable": {k: (None if e is None else str(sp.N(e, 4))) for k, e in found.items()}}
    if kept is not None:
        rec["pysr_kept"] = kept
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, f"{method}_{rcol}_{fold}.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    print(method, rcol, fold, rec["seconds"], "s", rec["readable"])


# --- scoring ---------------------------------------------------------------------------
def predict(expr, X):
    if expr is None:
        return np.full(len(X), np.nan)
    f = sp.lambdify([SYMS[v] for v in NAMES], expr, "numpy")
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        got = np.asarray(f(*X.T), dtype=complex)
    got = np.broadcast_to(got, (len(X),))
    bad = np.abs(got.imag) > 1e-9 * (1 + np.abs(got.real))
    return np.where(bad, np.nan, got.real)


def structure(expr):
    """Known-law check: is the answer c * gap / R (tau absent)? Also the power of each
    input when the answer is one product of powers."""
    if expr is None:
        return {"is_c_gap_over_R": False, "c": np.nan, "powers": None}
    ratio = sp.simplify(sp.powsimp(expr * SYMS["R"] / SYMS["gap"], force=True))
    ok = not ratio.free_symbols and ratio.is_number
    c = float(ratio) if ok else np.nan
    e = sp.powsimp(sp.expand_power_base(expr, force=True), force=True)
    coef, mono = e.as_coeff_Mul()
    pw = mono.as_powers_dict() if mono != 1 else {}
    powers = None
    if not isinstance(e, sp.Add) and all(b.is_Symbol and x.is_Number for b, x in pw.items()):
        powers = {str(b): float(x) for b, x in pw.items()}
    return {"is_c_gap_over_R": bool(ok), "c": c, "powers": powers}


def metrics(pred, y):
    good = np.isfinite(pred) & (pred > 0)
    rel = np.abs(pred - y) / y
    return {"median_abs_rel_err": float(np.nanmedian(np.where(np.isfinite(pred), rel, np.nan)))
            if np.any(np.isfinite(pred)) else np.nan,
            "rms_log_err": float(np.sqrt(np.mean(np.log(pred[good] / y[good]) ** 2)))
            if good.any() else np.nan,
            "frac_invalid": float(1 - good.mean())}


def known_c(X, y):
    """Best coefficient for the known form c*gap/R on these cells (least squares in log)."""
    return float(np.exp(np.mean(np.log(y * X[:, 0] / X[:, 2]))))


def score(out):
    import pandas as pd
    recs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(out, "*.json")))]
    rows = []
    for rcol in sorted({r["rcol"] for r in recs}):
        X, y = load(rcol)
        fl = folds(len(y))
        for fold in ["all"] + [str(i) for i in range(N_FOLDS)]:
            te = np.arange(len(y)) if fold == "all" else fl[int(fold)]
            tr = np.arange(len(y)) if fold == "all" else np.setdiff1d(np.arange(len(y)), te)
            c_fit = known_c(X[tr], y[tr])
            law_fit = c_fit * X[te, 2] / X[te, 0]
            base = [("known law, c = 1000 (physical)", PHYSICAL_C * X[te, 2] / X[te, 0], PHYSICAL_C),
                    ("known law, c fitted", law_fit, c_fit)]
            for name, pred, c in base:
                rows.append({"rcol": rcol, "fold": fold, "method": name,
                             "formula": f"{c:.4g}*gap/R", "is_c_gap_over_R": True, "c": c,
                             "powers": None, "seconds": 0.0,
                             "closeness_to_known_law": float(np.median(np.abs(pred - law_fit) / law_fit)),
                             **metrics(pred, y[te])})
            for r in recs:
                if r["rcol"] != rcol or r["fold"] != fold:
                    continue
                for pick, srep in r["formulas"].items():
                    expr = None if srep is None else sp.sympify(srep)
                    if expr is not None:
                        expr = expr.xreplace({s: SYMS[s.name] for s in expr.free_symbols})
                    pred = predict(expr, X[te])
                    st = structure(expr)
                    rows.append({"rcol": rcol, "fold": fold, "method": pick,
                                 "formula": r["readable"][pick], **st,
                                 "seconds": r["seconds"],
                                 "closeness_to_known_law": float(np.nanmedian(
                                     np.abs(pred - law_fit) / law_fit)),
                                 **metrics(pred, y[te])})
    df = pd.DataFrame(rows)
    df["c_vs_fitted"] = np.nan
    df["c_vs_physical"] = df["c"] / PHYSICAL_C
    for (rcol, fold), g in df.groupby(["rcol", "fold"]):
        cf = g.loc[g.method == "known law, c fitted", "c"].iloc[0]
        df.loc[g.index, "c_vs_fitted"] = g["c"] / cf
    df["exact_known_law"] = df.is_c_gap_over_R & (np.abs(df.c_vs_fitted - 1) <= COEF_TOL)
    df.to_csv(os.path.join(out, "allen_scores.csv"), index=False)
    return df


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    j = sub.add_parser("job")
    j.add_argument("--method", choices=list(HYBRID_VARIANTS) + ["PySINDy", "PySR"], required=True)
    j.add_argument("--rcol", choices=["ri", "input_resistance_mohm"], default="ri")
    j.add_argument("--fold", default="all")
    j.add_argument("--out", default=os.path.join(HERE, "results"))
    s = sub.add_parser("score")
    s.add_argument("--out", default=os.path.join(HERE, "results"))
    a = p.parse_args()
    if a.cmd == "job":
        job(a.method, a.rcol, a.fold, a.out)
    else:
        print(score(a.out).to_string())
