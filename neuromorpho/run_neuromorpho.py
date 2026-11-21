"""NeuroMorpho.org: neuron volume from surface area, total length and average diameter.
Hybrid vs PySINDy vs PySR, same method settings as allen/run_allen.py.

Data: NeuroMorpho.org API /api/morphometry (300,495 neurons), merged into
data/neuromorpho/morphometry.csv.
    S = surface (um^2), L = length (um), d = diameter (average, um);  V = volume (um^3)
Reference (cylinders): V = S*d/4 = (pi/4) d^2 L = S^2/(4 pi L). Equal for one cylinder;
on a whole neuron V >= S*d/4, because d is an average and V depends on d^2.

Splits (seed 0): one random permutation of all neurons. The first 50,000 are the held-out
test set; fits 0-4 each train on their own 2,000 neurons, disjoint, from the rest.
    python run_neuromorpho.py job --method Hybrid --fit 0
    python run_neuromorpho.py score
Every number is the method's formula evaluated on the held-out neurons' S, L, d and
compared with their recorded V.
"""
import argparse
import glob
import json
import os
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "allen"))

import sympy as sp  # noqa: E402

import run_allen as ra  # noqa: E402  (method code: run_hybrid, run_pysindy, run_pysr)

NAMES = ["S", "L", "d"]
SYMS = {v: sp.Symbol(v, positive=True) for v in NAMES}
N_TEST, N_TRAIN, N_FITS, SEED = 50_000, 2_000, 5, 0
OUT = os.path.join(HERE, "results")


def load():
    import pandas as pd
    df = pd.read_csv(os.path.join(REPO, "data", "neuromorpho", "morphometry.csv"))
    df = df.sort_values("neuron_id")
    X = df[["surface", "length", "diameter"]].to_numpy(float)
    y = df["volume"].to_numpy(float)
    ok = np.all(X > 0, axis=1) & (y > 0) & np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    return X[ok], y[ok]


def splits(n):
    perm = np.random.default_rng(SEED).permutation(n)
    test = perm[:N_TEST]
    trains = [perm[N_TEST + i * N_TRAIN: N_TEST + (i + 1) * N_TRAIN] for i in range(N_FITS)]
    return test, trains


def job(method, fit):
    X, y = load()
    tr = splits(len(y))[1][fit]
    t0 = time.time()
    kept = None
    if method.startswith("Hybrid"):
        found = ra.run_hybrid(X[tr], y[tr], names=NAMES, variant=method)
    elif method == "PySINDy":
        found = ra.run_pysindy(X[tr], y[tr], names=NAMES)
    else:
        found, kept = ra.run_pysr(X[tr], y[tr], names=NAMES)
    rec = {"method": method, "fit": fit, "n_train": int(len(tr)),
           "seconds": round(time.time() - t0, 2),
           "formulas": {k: (None if e is None else sp.srepr(e)) for k, e in found.items()},
           "readable": {k: (None if e is None else str(sp.N(e, 4))) for k, e in found.items()}}
    if kept is not None:
        rec["pysr_kept"] = kept
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"{method}_{fit}.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    print(method, fit, rec["seconds"], "s", rec["readable"])


def predict(expr, X):
    if expr is None:
        return np.full(len(X), np.nan)
    f = sp.lambdify([SYMS[v] for v in NAMES], expr, "numpy")
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        got = np.broadcast_to(np.asarray(f(*X.T), dtype=complex), (len(X),))
    bad = np.abs(got.imag) > 1e-9 * (1 + np.abs(got.real))
    return np.where(bad, np.nan, got.real)


REF_FORMS = {"c*S*d": SYMS["S"] * SYMS["d"], "c*d^2*L": SYMS["d"] ** 2 * SYMS["L"],
             "c*S^2/L": SYMS["S"] ** 2 / SYMS["L"]}


def form(expr):
    """Which cylinder form (if any) the answer is, and its coefficient."""
    if expr is None:
        return "", np.nan
    for name, base in REF_FORMS.items():
        r = sp.simplify(sp.powsimp(expr / base, force=True))
        if not r.free_symbols and r.is_number:
            return name, float(r)
    return "", np.nan


def size(expr):
    if expr is None:
        return np.nan, np.nan
    return len(expr.atoms(sp.Float)), int(sp.count_ops(expr))


def metrics(pred, y):
    good = np.isfinite(pred) & (pred > 0)
    rel = np.abs(pred - y) / y
    fin = np.isfinite(pred)
    return {"median_abs_rel_err": float(np.median(rel[fin])) if fin.any() else np.nan,
            "within_10pct": float(np.mean(np.where(fin, rel, np.inf) <= 0.10)),
            "rms_log_err": float(np.sqrt(np.mean(np.log(pred[good] / y[good]) ** 2)))
            if good.any() else np.nan,
            "frac_invalid": float(1 - good.mean())}


def score():
    import pandas as pd
    X, y = load()
    test, trains = splits(len(y))
    Xt, yt = X[test], y[test]
    rows = []
    ref = 0.25 * Xt[:, 0] * Xt[:, 2]
    for fit in range(N_FITS):
        tr = trains[fit]
        c_fit = float(np.exp(np.mean(np.log(y[tr] / (X[tr, 0] * X[tr, 2])))))
        for name, pred, c in [("reference S*d/4 (c = 0.25)", ref, 0.25),
                              ("reference c*S*d, c fitted", c_fit * Xt[:, 0] * Xt[:, 2], c_fit)]:
            rows.append({"fit": fit, "method": name, "formula": f"{c:.4g}*S*d",
                         "cylinder_form": "c*S*d", "c": c, "n_constants": 1, "n_ops": 2,
                         "seconds": 0.0, **metrics(pred, yt)})
    for p in sorted(glob.glob(os.path.join(OUT, "*_[0-9].json"))):
        r = json.load(open(p))
        for pick, srep in r["formulas"].items():
            expr = None if srep is None else sp.sympify(srep)
            if expr is not None:
                expr = expr.xreplace({s: SYMS[s.name] for s in expr.free_symbols})
            fm, c = form(expr)
            nc, no = size(expr)
            rows.append({"fit": r["fit"], "method": pick, "formula": r["readable"][pick],
                         "cylinder_form": fm, "c": c, "n_constants": nc, "n_ops": no,
                         "seconds": r["seconds"], **metrics(predict(expr, Xt), yt)})
    df = pd.DataFrame(rows).sort_values(["method", "fit"])
    df.to_csv(os.path.join(OUT, "neuromorpho_scores.csv"), index=False)
    return df


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    j = sub.add_parser("job")
    j.add_argument("--method", choices=list(ra.HYBRID_VARIANTS) + ["PySINDy", "PySR"], required=True)
    j.add_argument("--fit", type=int, required=True)
    sub.add_parser("score")
    a = p.parse_args()
    if a.cmd == "job":
        job(a.method, a.fit)
    else:
        print(score().drop(columns="formula").to_string())
