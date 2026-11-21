"""NeuroMorpho.org: total dendrite length from the bounding box and the number of branch
points (Cuntz, Mathy & Hausser 2012, PNAS: optimal wiring gives L ~ V^(1/3) n^(2/3)).
Hybrid vs PySINDy vs PySR, same method settings as allen/run_allen.py.

    W, H, D = width, height, depth of the neuron (um), given as separate inputs
    n       = n_bifs, number of branch points
    L       = length, total (um)
Neurons: depth >= 5 um (3D tracings) and n_bifs >= 1, from data/neuromorpho/morphometry.csv.
Reference: Cuntz's L = (3/(4 pi))^(1/3) V^(1/3) n^(2/3) with V = W*H*D (a bounding box, a
crude stand-in for the paper's spanning volume).

Splits (seed 0), as in run_neuromorpho.py: 50,000 held-out neurons; fits 0-4 on their own
2,000 neurons each.
    python run_wiring.py job --method Hybrid --fit 0
    python run_wiring.py score
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

import run_allen as ra  # noqa: E402

NAMES = ["W", "H", "D", "n"]
SYMS = {v: sp.Symbol(v, positive=True) for v in NAMES}
N_TEST, N_TRAIN, N_FITS, SEED = 50_000, 2_000, 5, 0
MIN_DEPTH = 5.0
OUT = os.path.join(HERE, "results_wiring")
CUNTZ_C = (3 / (4 * np.pi)) ** (1 / 3)


def load():
    import pandas as pd
    df = pd.read_csv(os.path.join(REPO, "data", "neuromorpho", "morphometry.csv"))
    df = df.sort_values("neuron_id")
    df = df[(df.depth >= MIN_DEPTH) & (df.width > 0) & (df.height > 0)
            & (df.n_bifs >= 1) & (df.length > 0)]
    X = df[["width", "height", "depth", "n_bifs"]].to_numpy(float)
    y = df["length"].to_numpy(float)
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    return X[ok], y[ok]


def splits(n):
    perm = np.random.default_rng(SEED).permutation(n)
    return perm[:N_TEST], [perm[N_TEST + i * N_TRAIN: N_TEST + (i + 1) * N_TRAIN]
                           for i in range(N_FITS)]


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


def powers(expr):
    """{input: exponent} if the answer is one product of powers, else None."""
    if expr is None:
        return None
    e = sp.powsimp(sp.expand_power_base(expr, force=True), force=True)
    if isinstance(e, sp.Add):
        return None
    _, mono = e.as_coeff_Mul()
    pw = mono.as_powers_dict() if mono != 1 else {}
    if not all(b.is_Symbol and x.is_Number for b, x in pw.items()):
        return None
    return {v: round(float(pw.get(SYMS[v], 0)), 4) for v in NAMES}


def metrics(pred, y):
    good = np.isfinite(pred) & (pred > 0)
    rel = np.abs(pred - y) / y
    fin = np.isfinite(pred)
    return {"median_abs_rel_err": float(np.median(rel[fin])) if fin.any() else np.nan,
            "within_10pct": float(np.mean(np.where(fin, rel, np.inf) <= 0.10)),
            "within_25pct": float(np.mean(np.where(fin, rel, np.inf) <= 0.25)),
            "rms_log_err": float(np.sqrt(np.mean(np.log(pred[good] / y[good]) ** 2)))
            if good.any() else np.nan,
            "frac_invalid": float(1 - good.mean())}


def score():
    import pandas as pd
    X, y = load()
    test, trains = splits(len(y))
    Xt, yt = X[test], y[test]
    cuntz_t = (Xt[:, 0] * Xt[:, 1] * Xt[:, 2]) ** (1 / 3) * Xt[:, 3] ** (2 / 3)
    rows = []
    for fit, tr in enumerate(trains):
        cuntz_tr = (X[tr, 0] * X[tr, 1] * X[tr, 2]) ** (1 / 3) * X[tr, 3] ** (2 / 3)
        c_fit = float(np.exp(np.mean(np.log(y[tr] / cuntz_tr))))
        A = np.column_stack([np.ones(len(tr)), np.log(X[tr])])
        k, *_ = np.linalg.lstsq(A, np.log(y[tr]), rcond=None)
        free = np.exp(k[0] + np.log(Xt) @ k[1:])
        free_txt = f"{np.exp(k[0]):.4g}*" + "*".join(f"{v}^{p:.3f}" for v, p in zip(NAMES, k[1:]))
        for name, pred, txt, pw in [
                ("reference Cuntz, textbook c", CUNTZ_C * cuntz_t,
                 f"{CUNTZ_C:.4g}*(W*H*D)^(1/3)*n^(2/3)", {"W": 1/3, "H": 1/3, "D": 1/3, "n": 2/3}),
                ("reference Cuntz, c fitted", c_fit * cuntz_t,
                 f"{c_fit:.4g}*(W*H*D)^(1/3)*n^(2/3)", {"W": 1/3, "H": 1/3, "D": 1/3, "n": 2/3}),
                ("reference free power law (log fit)", free, free_txt,
                 {v: round(float(p), 4) for v, p in zip(NAMES, k[1:])})]:
            rows.append({"fit": fit, "method": name, "formula": txt, "powers": pw,
                         "n_constants": 1 if "Cuntz" in name else 5, "seconds": 0.0,
                         **metrics(pred, yt)})
    for p in sorted(glob.glob(os.path.join(OUT, "*_[0-9].json"))):
        r = json.load(open(p))
        for pick, srep in r["formulas"].items():
            expr = None if srep is None else sp.sympify(srep)
            if expr is not None:
                expr = expr.xreplace({s: SYMS[s.name] for s in expr.free_symbols})
            rows.append({"fit": r["fit"], "method": pick, "formula": r["readable"][pick],
                         "powers": powers(expr),
                         "n_constants": np.nan if expr is None else len(expr.atoms(sp.Float)),
                         "seconds": r["seconds"], **metrics(predict(expr, Xt), yt)})
    df = pd.DataFrame(rows).sort_values(["method", "fit"])
    df.to_csv(os.path.join(OUT, "wiring_scores.csv"), index=False)
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
