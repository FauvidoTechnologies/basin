"""Recompute the headline numbers in busy_beaver/README.md from the committed files in results/.

    python busy_beaver/analyze.py            # all sections
    python busy_beaver/analyze.py classes    # one section: classes, speeds, bouncers, counters,
                                             #   champion, cryptids, halters

Needs numpy and scipy (installed with basin). Takes about a minute (most of it the counter file).
"""
import csv, gzip, json, math, os, re, sys, collections
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def load(name):
    with open(os.path.join(RES, name), newline="") as f:
        return list(csv.DictReader(f))


def classes():
    print("== Growth classes (window t = 2^12 .. 2^22) ==")
    order = ["1", "0", "1/2", "log", "1/3", "other"]
    for name in ("bb3_census.csv", "bb4_census.csv", "bb5_census.csv"):
        rows = load(name)
        for grp in sorted(set(r["group"] for r in rows)):
            sub = [r for r in rows if r["group"] == grp]
            c = collections.Counter(r["growth_class"] for r in sub)
            print(f"  {name[:3].upper()} {grp:6s} n={len(sub):6d}  " +
                  "  ".join(f"{k}: {100 * c[k] / len(sub):.2f}%" for k in order if c[k]))


def speeds():
    print("== Linear machines: speed of the growing edge (BB(5) sample) ==")
    v = [float(r["rate"]) for r in load("bb5_census.csv") if r["growth_class"] == "1"]
    def frac(x):
        f = Fraction(x).limit_denominator(12)
        return f if abs(float(f) - x) < 2e-3 else None
    c = collections.Counter(str(frac(x)) if frac(x) is not None else "irregular" for x in v)
    print("  most common:", c.most_common(14))
    # parity rule: d cells in p steps needs p = d (mod 2); in lowest terms a/b with a+b odd needs (2a, 2b)
    even = sum(n for k, n in c.items() if k not in ("irregular", "None") and (Fraction(k).numerator + Fraction(k).denominator) % 2 == 0)
    odd = sum(n for k, n in c.items() if k not in ("irregular", "None") and (Fraction(k).numerator + Fraction(k).denominator) % 2 == 1)
    print(f"  speeds a/b with a+b even (period can be b): {even};  a+b odd (period must be 2b): {odd}")


def bouncers():
    print("== sqrt(t) machines: is c^2 a simple fraction? (BB(5) sample) ==")
    c2 = np.array([float(r["rate"]) ** 2 for r in load("bb5_census.csv") if r["growth_class"] == "1/2"])
    rng = np.random.default_rng(1)
    null = np.exp(rng.uniform(np.log(c2.min()), np.log(c2.max()), 20000))
    hit = lambda xs, den, tol: np.mean([abs(float(Fraction(x).limit_denominator(den)) - x) <= tol * x for x in xs])
    for den, tol in [(12, .004), (6, .002), (4, .001)]:
        print(f"  denominator <= {den:2d}, within {100 * tol:.1f}%:  machines {100 * hit(c2, den, tol):5.1f}%   random {100 * hit(null, den, tol):5.1f}%")
    common = collections.Counter(str(Fraction(x).limit_denominator(6)) for x in c2 if abs(float(Fraction(x).limit_denominator(6)) - x) <= .002 * x)
    print("  most common c^2:", common.most_common(10))
    for name in ("bb4_census.csv", "bb5_census.csv"):
        third = [r for r in load(name) if r["growth_class"] == "1/3"]
        print(f"  {name[:3].upper()} t^(1/3) machines: " + ", ".join(f"{r['machine']} c^3={float(r['rate']) ** 3:.2f}" for r in third))


def counters():
    print("== Counter bases: all 297,225 counter candidates in the BB(5) time-limit block ==")
    with gzip.open(os.path.join(RES, "bb5_counter_bases.csv.gz"), "rt", newline="") as f:
        rows = list(csv.DictReader(f))
    c = collections.Counter(r["base"] or "too few increments" for r in rows)
    fp = collections.defaultdict(set)
    for r in rows: fp[r["base"] or "too few increments"].add(r["fingerprint"])
    for k, n in c.most_common():
        print(f"  {k:20s} {n:8,d} records  {100 * n / len(rows):6.2f}%   {len(fp[k]):7,d} distinct growth sequences")
    F = np.array([float(r["time_factor_per_new_cell"]) for r in rows])
    rr = lambda co: max(z.real for z in np.roots(co) if abs(z.imag) < 1e-12)
    phi = (1 + 5 ** .5) / 2
    ref = {"phi": phi, "plastic": rr([1, 0, -1, -1]), "plastic^2": rr([1, 0, -1, -1]) ** 2, "tribonacci": rr([1, -1, -1, -1]),
           "narayana": rr([1, -1, 0, -1]), "sqrt(1+sqrt2)": (1 + 2 ** .5) ** .5}
    for k, v in ref.items():
        e = np.abs(F / v - 1); m = e < 1e-3
        print(f"  cluster at {k:14s} = {v:.6f}: {int(m.sum()):6,d} records, median relative error {np.median(e[m]):.1e}")
    un = collections.Counter(np.round(F[[r["base"] == "unidentified" for r in rows]], 3))
    print("  largest unidentified clusters:", un.most_common(6))


def champion():
    print("== BB(5) champion: span/sqrt(t) from its Collatz-like rule vs simulation ==")
    lines = open(os.path.join(RES, "champion_trace.txt")).read().splitlines()
    cfg = [(int(l.split()[1][2:]), int(l.split()[5][5:])) for l in lines if l.startswith("CFG")]
    print("  all 14 configurations have the form 0^inf <A 1^n 0^inf:", all("form_ok=1" in l for l in lines if l.startswith("CFG")))
    print(f"  predicted at configurations 4/sqrt(5) = {4 / math.sqrt(5):.4f}; measured:",
          [round(s / math.sqrt(t), 4) for t, s in cfg[8:]])
    S = np.array([[int(x) for x in l.split()[1:]] for l in lines if l.startswith("S ")], float)
    t, r = S[:, 0], S[:, 1] / np.sqrt(S[:, 0])
    for (t0, _), (t1, _) in zip(cfg[8:12], cfg[9:13]):
        m = (t >= t0) & (t < t1); i = np.argmax(r[m])
        print(f"  cycle x{t1 / t0:.3f}: min {r[m].min():.4f}, peak {r[m][i]:.4f} at {100 * math.log(t[m][i] / t0) / math.log(t1 / t0):.0f}% (log t)")
    tau = np.linspace(0, 5 / 9, 2_000_001); tt = 5 / 16 + tau; R = (1 + 0.4 * np.sqrt(5 * tau)) / np.sqrt(tt)
    w = np.gradient(np.log(tt))
    print(f"  model: min {R.min():.4f} (4/sqrt5), max {R.max():.4f} (2), mean over log t {(R * w).sum() / w.sum():.4f}; basin fit c = 1.920")


def cryptids():
    from scipy.stats import mannwhitneyu
    print("== BB(6): roughness (cross-validated error) within growth class ==")
    rows = load("bb6_fits.csv"); bb5 = load("bb5_census.csv")
    for r in rows: r["cv"] = float(r["cv_rel_mse"])
    ref = collections.defaultdict(list)
    for r in rows:
        if r["set"] in ("holdout", "random6"): ref[(r["set"], r["growth_class"])].append(r["cv"])
    for r in bb5: ref[("bb5", r["growth_class"])].append(float(r["cv_rel_mse"]))
    pct = lambda g, k, v: (100 * np.mean(np.array(ref[(g, k)]) < v)) if len(ref[(g, k)]) >= 10 else float("nan")
    for g in ("random6", "holdout"):
        c = collections.Counter(r["growth_class"] for r in rows if r["set"] == g); n = sum(c.values())
        print(f"  {g:8s} n={n}: " + ", ".join(f"{k} {100 * v / n:.1f}%" for k, v in c.most_common()))
    nb = lambda xs: np.array([x for x in xs if x > 0])
    ho = nb([r["cv"] for r in rows if r["set"] == "holdout" and r["growth_class"] != "0"])
    for name, refv in (("random 6-state", nb([r["cv"] for r in rows if r["set"] == "random6" and r["growth_class"] != "0"])),
                       ("BB(5) census", nb([float(r["cv_rel_mse"]) for r in bb5 if r["growth_class"] != "0"]))):
        auc = mannwhitneyu(ho, refv).statistic / (len(ho) * len(refv))
        print(f"  holdouts rougher than the {name} 99th percentile: {100 * np.mean(ho > np.percentile(refv, 99)):.0f}%  (AUC {auc:.3f})")
    print(f"  {'machine':18s} {'status':26s} class  {'cv':>8s}  pct-in-holdouts pct-in-random6 pct-in-BB(5)")
    for r in rows:
        if r["set"] in ("holdout", "random6"): continue
        k = r["growth_class"]
        print(f"  {r['name'][:18]:18s} {r['status'][:26]:26s} {k:5s} {r['cv']:8.1e}  {pct('holdout', k, r['cv']):6.1f}        "
              f"{pct('random6', k, r['cv']):6.1f}        {pct('bb5', k, r['cv']):6.1f}")
    cr = [r for r in rows if r["set"] == "cryptid"]
    H = [pct("holdout", r["growth_class"], r["cv"]) for r in cr if r["status"] == "probviously halting"]
    N = [pct("holdout", r["growth_class"], r["cv"]) for r in cr if r["status"] == "probviously non-halting"]
    print(f"  halting cryptids median {np.median(H):.0f} (n={len(H)}), non-halting median {np.median(N):.0f} (n={len(N)}),"
          f" Mann-Whitney p = {mannwhitneyu(H, N).pvalue:.2f}")
    FL = re.compile(r"(?<![\w.])-?\d+\.\d+(e[-+]?\d+)?")
    clean = {"1/2": {"c*sqrt(t)", "c*sqrt(t) + c"}, "1/3": {"c*t**(1/3)", "c*t**(1/3) + c"}, "log": {"c*u", "c*u + c"}}
    for k, ok in clean.items():
        xs = [r for r in rows if r["set"] == "holdout" and r["growth_class"] == k]
        n = sum(FL.sub("c", r["law"]).replace(" - ", " + ").replace("c*(", "(") in ok for r in xs)
        print(f"  holdouts in class {k}: {len(xs)}, with a clean law {n} ({100 * n / max(len(xs), 1):.0f}%)")


def halters():
    from scipy.stats import mannwhitneyu
    import sympy as sp
    print("== BB(5) halters running 2^12..2^16 steps vs non-halters, same window 2^10..2^12 ==")
    T = sp.Symbol("t", positive=True); FL = re.compile(r"(?<![\w.])-?\d+\.\d+(e[-+]?\d+)?"); nco = lambda l: len(FL.findall(l))
    NICE = [(0, "0"), (1 / 3, "1/3"), (1 / 2, "1/2"), (2 / 3, "2/3"), (1, "1")]
    out = collections.defaultdict(list)
    for line in open(os.path.join(RES, "halters_vs_nonhalters_window_2e10_2e12.jsonl")):
        d = json.loads(line)
        if d.get("cls"): out[d["set"]].append((d["cls"], 0.0)); continue
        lt, lu, ct, cu = d["t"]["law"], d["u"]["law"], d["t"]["cv"], d["u"]["cv"]
        use_u = bool(lu) and ((nco(lu) < nco(lt)) if max(ct, cu) <= 2 * min(ct, cu) and nco(lu) != nco(lt) else cu < ct)
        if use_u: out[d["set"]].append(("log", cu)); continue
        f = sp.lambdify(T, sp.sympify(lt, locals={"t": T}), "numpy"); x = 2.0 ** 12
        with np.errstate(all="ignore"): a = float(np.log(f(x * 1.1) / f(x / 1.1)) / np.log(1.21))
        near = min(NICE, key=lambda n: abs(n[0] - a)) if np.isfinite(a) else (None, "?")
        out[d["set"]].append((near[1] if np.isfinite(a) and abs(near[0] - a) < 0.04 else "other", ct))
    for s in ("halter", "nonhalt"):
        c = collections.Counter(k for k, _ in out[s]); n = len(out[s])
        print(f"  {s:8s} n={n:5d} " + ", ".join(f"{k} {100 * v / n:.1f}%" for k, v in c.most_common()))
    for k in ("1/2", "other"):
        h = [v for c, v in out["halter"] if c == k]; nh = [v for c, v in out["nonhalt"] if c == k]
        print(f"  class {k}: median cv halters {np.median(h):.1e} (n={len(h)}) vs non-halters {np.median(nh):.1e} (n={len(nh)}),"
              f" Mann-Whitney p = {mannwhitneyu(h, nh).pvalue:.1g}")


SECTIONS = dict(classes=classes, speeds=speeds, bouncers=bouncers, counters=counters,
                champion=champion, cryptids=cryptids, halters=halters)
if __name__ == "__main__":
    for name in (sys.argv[1:] or SECTIONS):
        SECTIONS[name](); print()
