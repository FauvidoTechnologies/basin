"""Hybrid law finder: BACON proposes terms, sparse regression decides (Test 13 on).

Built from what Tests 1-12 showed: PySINDy decides well (one fit on all the data) but can
only find laws inside its fixed term list; BACON can build the terms real laws need
(powers like 2/3, sums like L + 3 in a denominator) but decides with yes/no tolerance
tests on small slices, which break under noise.

    target y, inputs x_1..x_n (all rows used at every step)

1. Propose powers, BACON's way. For each input x_j, use the experiment's structure: rows
   that share the values of every other input form a slice. Fit
       log|y| = (one intercept per slice) + p_j log x_j
   on all slices at once, which gives p_j and its standard error. If p_j is within
   max(2 SE, 0.02) of a simple fraction (denominator 1, 2 or 3, |p| <= 3, not 0), x_j^p
   is proposed (3 SE, floor 0.03; the 2 nearest such fractions at most). Nothing is accepted or rejected here; proposals only add candidate terms.
   (If the data has no slices, as in a time series, one fit of log|y| on all log x_j.)
2. Term list: every product of non-negative integer powers with total degree <= 2
   (1, x_i, x_i x_j, x_i^2), plus each proposed x_j^p, plus the product of all proposed
   powers (BACON's invariant, e.g. T*D^2).
3. Implicit form, as in SINDy-PI: y = N(x)/D(x), with N and D sparse sums of terms, fitted
   as the linear problem  y*D(x) = N(x),  with one D term's coefficient set to 1 (each
   term in turn). Ratio laws like (M1 T1 + M2 T2)/(M1 + M2) or 3/(X + 2) become sums.
   Candidate laws: every law with up to 4 terms in total (screened by a quick fit on all
   rows; the best 30 of each size go on), plus larger ones from two searches per D term:
   PySINDy's STLSQ over the usual threshold grid (normalized columns, no ridge penalty),
   and forward selection (add the term that most reduces the squared error, up to
   MAX_FREE free coefficients).
4. Decide by held-out error, not tolerances. Every candidate with at most MAX_FREE free
   coefficients is scored by 5-fold cross-validation repeated over 5 random splits
   (leave-one-out below 10 rows): fit on
   the other folds, predict y = N/D on the held-out rows, mean squared relative error.
   Coefficients are first fitted linearly, then refined by least squares on the
   relative prediction error itself. Pick the fewest terms whose error is within one
   standard error of the best (the "one-standard-error rule"; differences below 1e-10
   count as rounding), ties going to lower total power.
Not built (yet): a further round that proposes new terms from what the chosen law leaves
unexplained.
"""
from fractions import Fraction
from itertools import combinations_with_replacement

import numpy as np
import pysindy as ps
import sympy as sp
import warnings

from scipy.optimize import least_squares

# PySINDy STLSQ threshold grid on normalized columns, copied from bacon-vs-pysindy's
# sindy_runner.THRESHOLDS_NORM (same values) so this file stands alone.
THRESHOLDS_NORM = np.logspace(-3, 3, 25)

MAX_FREE = 5            # most free coefficients a law may have (N and D together)
EXHAUSTIVE_TERMS = 4    # every law with up to this many terms (N and D together) is tried
SCREEN_KEEP = 30        # per law size, how many survive the quick screen into cross-validation
POWER_SE = 3            # step 1: accept a fraction within POWER_SE standard errors ...
POWER_FLOOR = 0.03      # ... or within this, whichever is larger (clean data has SE ~0)
POWER_KEEP = 2          # ... at most the 2 nearest fractions per input
CV_FLOOR = 1e-10        # errors below this (relative RMS 1e-5) count as equal: rounding
FOLDS = 5
REPEATS = 5             # repeated cross-validation: 5 different random splits into folds
CV_SEED = 0
RATIONALS = sorted({Fraction(k, d) for d in (1, 2, 3) for k in range(-3 * d, 3 * d + 1)
                    if k != 0})
PARAMS = {"max_free": MAX_FREE, "power_se": POWER_SE, "power_floor": POWER_FLOOR,
          "cv_floor": CV_FLOOR, "folds": FOLDS, "repeats": REPEATS, "exhaustive_terms": EXHAUSTIVE_TERMS,
          "screen_keep": SCREEN_KEEP, "refine": "least squares on relative error",
          "cv_seed": CV_SEED, "power_keep": POWER_KEEP, "power_denominators": [1, 2, 3], "max_abs_power": 3,
          "integer_degree": 2, "stlsq_grid": "THRESHOLDS_NORM (25, 0.001-1000)",
          "selection": "one-standard-error rule on cross-validated relative error"}


# --- 1. propose powers --------------------------------------------------------------
def propose_powers(X, y):
    """{input index: Fraction power} for inputs whose power law fits a simple fraction."""
    n_rows, n = X.shape
    if np.any(y == 0):
        return {}, {}
    ly = np.log(np.abs(y))
    positive = [j for j in range(n) if np.all(X[:, j] > 0)]
    estimates, out = {}, {}
    for j in positive:
        others = np.delete(X, j, axis=1)
        if others.shape[1]:
            _, sid = np.unique(others, axis=0, return_inverse=True)
        else:
            sid = np.zeros(n_rows, dtype=int)
        k = sid.max() + 1
        if n_rows - k - 1 >= 1 and k <= n_rows / 2:      # slices exist: BACON's design
            A = np.column_stack([np.eye(k)[sid], np.log(X[:, j])])
            col = A.shape[1] - 1
        else:                                            # no slices: one global fit
            A = np.column_stack([np.ones(n_rows)] + [np.log(X[:, i]) for i in positive])
            col = 1 + positive.index(j)
        coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
        dof = n_rows - np.linalg.matrix_rank(A)
        if dof < 1:
            continue
        resid = ly - A @ coef
        try:
            cov = (resid @ resid / dof) * np.linalg.pinv(A.T @ A)
        except np.linalg.LinAlgError:
            continue
        p, se = float(coef[col]), float(np.sqrt(max(cov[col, col], 0.0)))
        estimates[j] = (p, se)
        tol = max(POWER_SE * se, POWER_FLOOR)
        near = sorted((r for r in RATIONALS if abs(p - float(r)) <= tol),
                      key=lambda r: abs(p - float(r)))[:POWER_KEEP]
        if near:
            out[j] = near          # nearest first
    return out, estimates


# --- 2. term list -------------------------------------------------------------------
def term_list(n, powers):
    """Exponent vectors (tuples of Fractions), constant first, no duplicates."""
    terms = [tuple(Fraction(0) for _ in range(n))]
    for deg in (1, 2):
        for combo in combinations_with_replacement(range(n), deg):
            e = [Fraction(0)] * n
            for i in combo:
                e[i] += 1
            terms.append(tuple(e))
    for j, ps_ in powers.items():
        for p in ps_:
            e = [Fraction(0)] * n
            e[j] = p
            terms.append(tuple(e))
    if len(powers) >= 2:          # BACON's invariant, from each input's nearest fraction
        e = [Fraction(0)] * n
        for j, ps_ in powers.items():
            e[j] = ps_[0]
        terms.append(tuple(e))
    seen, out = set(), []
    for t in terms:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def columns(X, terms):
    cols = []
    for e in terms:
        c = np.ones(len(X))
        for j, p in enumerate(e):
            if p != 0:
                c = c * (X[:, j] ** int(p) if p.denominator == 1 else X[:, j] ** float(p))
        cols.append(c)
    return np.column_stack(cols)


# --- 3. implicit fit ----------------------------------------------------------------
def _design(theta, y, anchor):
    """Features for y*theta[anchor] = sum a_k theta_k - sum_{k != anchor} b_k y theta_k."""
    m = theta.shape[1]
    feats = [theta[:, k] for k in range(m)] + [-y * theta[:, k] for k in range(m) if k != anchor]
    keys = [("N", k) for k in range(m)] + [("D", k) for k in range(m) if k != anchor]
    return np.column_stack(feats), keys, y * theta[:, anchor]


def _lstsq(F, t):
    w, *_ = np.linalg.lstsq(F, t, rcond=None)
    return w


def _candidates(theta, y):
    m = theta.shape[1]
    found = {}
    for anchor in range(m):
        F, keys, t = _design(theta, y, anchor)
        supports = []
        # STLSQ path (PySINDy)
        for thr in THRESHOLDS_NORM:
            opt = ps.STLSQ(threshold=thr, normalize_columns=True, alpha=0.0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    opt.fit(F, t.reshape(-1, 1))
                except Exception:
                    continue
            s = tuple(np.flatnonzero(np.asarray(opt.coef_).ravel()))
            if 0 < len(s) <= MAX_FREE:
                supports.append(s)
        # forward selection
        chosen = []
        scale = np.linalg.norm(F, axis=0)
        scale[scale == 0] = 1
        Fn = F / scale
        for _ in range(MAX_FREE):
            best, best_rss = None, np.inf
            for k in range(F.shape[1]):
                if k in chosen:
                    continue
                cols = chosen + [k]
                w = _lstsq(Fn[:, cols], t)
                rss = float(np.sum((t - Fn[:, cols] @ w) ** 2))
                if rss < best_rss:
                    best, best_rss = k, rss
            if best is None:
                break
            chosen.append(best)
            supports.append(tuple(sorted(chosen)))
        for s in supports:
            num = frozenset(keys[i][1] for i in s if keys[i][0] == "N")
            den = frozenset([anchor] + [keys[i][1] for i in s if keys[i][0] == "D"])
            if num and (num, den) not in found:
                found[(num, den)] = (anchor, s)
    return found


def _unpack(w, keys, support, m, anchor):
    a, b = np.zeros(m), np.zeros(m)
    b[anchor] = 1.0
    for wi, i in zip(w, support):
        side, k = keys[i]
        (a if side == "N" else b)[k] = wi
    return a, b


def _exhaustive(theta, y):
    """Every (N, D) with |N| + |D| <= EXHAUSTIVE_TERMS, screened by the training error of
    the linear implicit fit; the best SCREEN_KEEP of each size are kept."""
    from itertools import combinations
    m = theta.shape[1]
    scale = np.abs(y) if (np.all(y > 0) or np.all(y < 0)) else np.std(y)
    by_size = {}
    for size in range(2, EXHAUSTIVE_TERMS + 1):
        for nd in range(1, size):
            for den in combinations(range(m), size - nd):
                anchor, rest = den[0], den[1:]
                t = y * theta[:, anchor]
                Dpart = -y[:, None] * theta[:, list(rest)]
                for num in combinations(range(m), nd):
                    F = np.column_stack([theta[:, list(num)], Dpart])
                    w, *_ = np.linalg.lstsq(F, t, rcond=None)
                    Nv = theta[:, list(num)] @ w[:nd]
                    Dv = theta[:, anchor] + (theta[:, list(rest)] @ w[nd:] if rest else 0)
                    with np.errstate(all="ignore"):
                        err = np.mean(((Nv / Dv - y) / scale) ** 2)
                    if np.isfinite(err):
                        by_size.setdefault(size, []).append((err, frozenset(num), frozenset(den)))
    out = {}
    for size, lst in by_size.items():
        lst.sort(key=lambda r: r[0])
        for _, num, den in lst[:SCREEN_KEEP]:
            out[(num, den)] = None
    return out


def _support_for(num, den, m):
    """(anchor, support indices into _design's features) for a given N and D."""
    anchor = min(den)
    keys = [("N", k) for k in range(m)] + [("D", k) for k in range(m) if k != anchor]
    sup = tuple(sorted([keys.index(("N", k)) for k in num] +
                       [keys.index(("D", k)) for k in den if k != anchor]))
    return anchor, sup


def _fit_predict(theta_tr, y_tr, theta_te, anchor, support):
    """Linear implicit fit, then refined on the prediction error itself: the implicit
    residual y*D - N weights rows by D and puts the noise on both sides, which biases the
    coefficients under noise."""
    F, keys, t = _design(theta_tr, y_tr, anchor)
    m = theta_tr.shape[1]
    w0 = _lstsq(F[:, list(support)], t)
    scale = np.abs(y_tr) if (np.all(y_tr > 0) or np.all(y_tr < 0)) else np.std(y_tr)

    def resid(w):
        a, b = _unpack(w, keys, support, m, anchor)
        with np.errstate(all="ignore"):
            r = ((theta_tr @ a) / (theta_tr @ b) - y_tr) / scale
        return np.where(np.isfinite(r), r, 1e6)
    try:
        w = least_squares(resid, w0, method="lm" if len(y_tr) >= len(w0) else "trf",
                          max_nfev=200 * (len(w0) + 1)).x
        if np.sum(resid(w) ** 2) > np.sum(resid(w0) ** 2):
            w = w0
    except Exception:
        w = w0
    a, b = _unpack(w, keys, support, m, anchor)
    with np.errstate(all="ignore"):
        pred = (theta_te @ a) / (theta_te @ b)
    return pred, a, b


def _errors(pred, y):
    if np.all(y > 0) or np.all(y < 0):
        e = (pred - y) / y
    else:
        e = (pred - y) / np.std(y)
    return float(np.mean(e ** 2)) if np.all(np.isfinite(e)) else np.inf


def _folds(n_rows):
    """REPEATS x FOLDS held-out sets (leave-one-out below 10 rows). Errors are averaged
    over all of them; the standard error uses one value per repeat."""
    if n_rows < 10:
        return [[np.array([i]) for i in range(n_rows)]]
    rng = np.random.default_rng(CV_SEED)
    reps = []
    for _ in range(REPEATS):
        idx = rng.permutation(n_rows)
        reps.append([idx[i::FOLDS] for i in range(FOLDS)])
    return reps


# --- 4. decide ----------------------------------------------------------------------
def discover(X, y, symbols, details=False):
    """Returns a sympy expression for y in terms of `symbols` (positive inputs assumed
    only where fractional powers are proposed)."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    powers, estimates = propose_powers(X, y)
    terms = term_list(X.shape[1], powers)
    theta = columns(X, terms)
    cands = _candidates(theta, y)
    for key in _exhaustive(theta, y):
        if key not in cands:
            cands[key] = _support_for(*key, theta.shape[1])
    folds = _folds(len(y))
    scored = []
    for (num, den), (anchor, support) in cands.items():
        per_rep = []
        for rep in folds:
            errs = []
            for te in rep:
                tr = np.setdiff1d(np.arange(len(y)), te)
                pred, _, _ = _fit_predict(theta[tr], y[tr], theta[te], anchor, support)
                errs.append(_errors(pred, y[te]))
            per_rep.append(errs)
        flat = np.array([e for r in per_rep for e in r])
        cv = float(np.mean(flat))
        if not np.all(np.isfinite(flat)):
            se = np.inf
        elif len(per_rep) > 1:   # repeats share rows, so SE from fold scores within a repeat
            se = float(np.mean([np.std(r, ddof=1) / np.sqrt(len(r)) for r in per_rep]))
        else:
            se = float(np.std(flat, ddof=1) / np.sqrt(len(flat)))
        degree = sum(float(sum(abs(p) for p in terms[k])) for k in num | den)
        scored.append((cv, se, len(num) + len(den), anchor, support, num, den, degree))
    finite = [s for s in scored if np.isfinite(s[0])]
    if not finite:
        return (None, {}) if details else None
    best = min(finite, key=lambda s: s[0])
    ok = [s for s in finite if s[0] <= best[0] + max(best[1], CV_FLOOR)]
    # fewest terms; then lowest total power (so u/(w+5), not u*w/(w^2+5w)); then error
    cv, se, nterms, anchor, support, num, den, _ = min(ok, key=lambda s: (s[2], s[7], s[0]))
    _, a, b = _fit_predict(theta, y, theta, anchor, support)

    def expr(e):
        return sp.Mul(*[s ** sp.Rational(p.numerator, p.denominator)
                        for s, p in zip(symbols, e) if p != 0])
    N = sum((float(c) * expr(e) for c, e in zip(a, terms) if c != 0), sp.Integer(0))
    D = sum((float(c) * expr(e) for c, e in zip(b, terms) if c != 0), sp.Integer(0))
    law = N / D
    if details:
        return law, {"proposed_powers": {str(symbols[j]): [str(p) for p in ps_] for j, ps_ in powers.items()},
                     "power_estimates": {str(symbols[j]): (round(p, 4), round(s, 4))
                                         for j, (p, s) in estimates.items()},
                     "n_terms_in_list": len(terms), "n_candidates": len(scored),
                     "cv_rel_mse": cv, "cv_se": se}
    return law
