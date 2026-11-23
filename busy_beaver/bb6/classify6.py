# Classify bb6_fits.jsonl with the same rules as census/classify.py.
import csv, json, re, collections, numpy as np, sympy as sp
T = sp.Symbol("t", positive=True); U = sp.Symbol("u", positive=True)
FLOAT = re.compile(r"(?<![\w.])-?\d+\.\d+(e[-+]?\d+)?")
ncoef = lambda law: len(FLOAT.findall(law))
NICE = [(0, "0"), (1/3, "1/3"), (1/2, "1/2"), (2/3, "2/3"), (1, "1")]
special = list(csv.DictReader(open("special.csv")))
names = {"holdout": [l.strip() for l in open("holdouts_m.txt")], "random6": [l.strip() for l in open("rand6.txt")]}
out = []
for line in open("bb6_fits.jsonl"):
    d = json.loads(line)
    if d["set"] == "special":
        s = special[d["i"]]; d.update(group=s["set"], name=s["name"], status=s["status"], machine=s["machine"])
    else:
        d.update(group=d["set"], name="", status="", machine=names[d["set"]][d["i"]])
    if d.get("cls") == "bounded":
        d.update(exp="0", law=str(d["span"]), cv=0.0); out.append(d); continue
    lt, lu, ct, cu = d["t"]["law"], d["u"]["law"], d["t"]["cv"], d["u"]["cv"]
    if lu and lt:
        close = max(ct, cu) <= 2 * min(ct, cu)
        use_u = (ncoef(lu) < ncoef(lt)) if close and ncoef(lu) != ncoef(lt) else cu < ct
    else:
        use_u = bool(lu)
    if use_u:
        d.update(exp="log", law=lu, cv=cu)
    else:
        e = sp.sympify(lt, locals={"t": T}); f = sp.lambdify(T, e, "numpy"); x = float(d["last_t"])
        with np.errstate(all="ignore"):
            a = float(np.log(f(x * 1.1) / f(x / 1.1)) / np.log(1.21))
        near = min(NICE, key=lambda n: abs(n[0] - a)) if np.isfinite(a) else (None, "?")
        d.update(exp=near[1] if np.isfinite(a) and abs(near[0] - a) < 0.04 else "other", law=lt, cv=ct, alpha=a)
    out.append(d)
json.dump(out, open("bb6_classified.json", "w"))
for g in ("random6", "holdout"):
    c = collections.Counter(d["exp"] for d in out if d["group"] == g); n = sum(c.values())
    print(f"{g:8s} n={n:5d}", {k: f"{v} ({100*v/n:.1f}%)" for k, v in c.most_common()})
