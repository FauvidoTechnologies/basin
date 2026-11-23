import sys, json, re, collections, numpy as np, sympy as sp
from sympy.parsing.sympy_parser import parse_expr
K = 49
SIM, FITS, BIN, IDX, OUT, NST = (sys.argv[1:7] if len(sys.argv) > 6 else
    ("sim22.out", "fits.jsonl", "sample.bin", "sample_idx.npy", "classified.json", "5"))
NST = int(NST)
r = np.fromfile(SIM, np.int32).reshape(-1, 1 + 2 * K)
idx = np.load(IDX); raw = open(BIN, "rb").read()
_sp = r[:, 1:1 + K].astype(float)
RAWSLOPE = np.log(_sp[:, -1] / np.maximum(_sp[:, -9], 1)) / np.log(4)  # 2^20..2^22
T = sp.Symbol("t", positive=True); U = sp.Symbol("u", positive=True)
def machine(i):
    b = raw[30*i:30*i+30]; out = []
    for s in range(NST):
        p = ""
        for q in range(2):
            w, m, g = b[6*s+3*q:6*s+3*q+3]
            p += "---" if g == 0 else f"{w}{'RL'[m]}{chr(64+g)}"
        out.append(p)
    return "_".join(out)
FLOAT = re.compile(r"(?<![\w.])-?\d+\.\d+(e[-+]?\d+)?")
def shape(law):
    k = iter(range(99))
    s = FLOAT.sub(lambda m: f"c{next(k)}", law)
    syms = {f"c{j}": sp.Symbol(f"c{j}") for j in range(99)}
    syms.update(t=T, u=U)
    try:
        s = str(sp.sympify(s, locals=syms))
    except Exception:
        pass
    s = re.sub(r"c\d+", "c", s).replace(" - ", " + ")
    s = re.sub(r"^-", "", s)
    s = re.sub(r"^c\*\((.*)\)(/[^()]+)?$", lambda m: "(" + m.group(1) + ")" + (m.group(2) or ""), s)
    s = re.sub(r"\(c \+ (c\*[^()+]+)\)", r"(\1 + c)", s)
    return s
def ncoef(law):
    return len(FLOAT.findall(law))
def slope(expr, sym, x):
    f = sp.lambdify(sym, expr, "numpy")
    with np.errstate(all="ignore"):
        a, b = float(f(x / 1.1)), float(f(x * 1.1))
    return np.log(b / a) / np.log(1.21) if a > 0 and b > 0 else np.nan
NICE = [(0, "0"), (1/3, "1/3"), (1/2, "1/2"), (2/3, "2/3"), (1, "1")]
rows = []
for line in open(FITS):
    d = json.loads(line); i = d["i"]
    d["machine"] = machine(i); d["group"] = ("time" if idx[i] < 14322029 else "space") if NST == 5 else "all"
    if d.get("cls") == "bounded":
        d.update(exp="0", law=str(d["span"]), form="bounded", rate=float(d["span"])); rows.append(d); continue
    ct, cu = d["t"]["cv"], d["u"]["cv"]
    lt, lu = d["t"]["law"], d["u"]["law"]
    if RAWSLOPE[i] > 0.4:
        lu = None  # a log law cannot describe growth this fast; same rule as small/fit_small.py
    if lu and lt:
        close = max(ct, cu) <= 2 * min(ct, cu)
        use_u = (ncoef(lu) < ncoef(lt)) if close and ncoef(lu) != ncoef(lt) else cu < ct
    else:
        use_u = bool(lu)
    if use_u:
        e = parse_expr(d["u"]["law"], {"u": U}); d.update(cls="log", law=d["u"]["law"], form=shape(d["u"]["law"]), cv=cu)
        d["exp"] = "log"
        f = sp.lambdify(U, e, "numpy"); d["rate"] = float(f(22.5) - f(21.5))  # cells per doubling of t
    else:
        e = parse_expr(d["t"]["law"], {"t": T}); d.update(law=d["t"]["law"], form=shape(d["t"]["law"]), cv=ct)
        a = slope(e, T, 2.0**22); d["alpha"] = a
        near = min(NICE, key=lambda n: abs(n[0] - a)) if np.isfinite(a) else (None, "?")
        d["exp"] = near[1] if np.isfinite(a) and abs(near[0] - a) < 0.04 else "other"
        d["cls"] = "power"
        if d["exp"] in ("1/3", "1/2", "2/3", "1"):
            p = {"1/3": 1/3, "1/2": 1/2, "2/3": 2/3, "1": 1.0}[d["exp"]]
            d["rate"] = float(sp.lambdify(T, e, "numpy")(2.0**22)) / 2.0**(22 * p)
    rows.append(d)
json.dump(rows, open(OUT, "w"))
N = len(rows)
print("machines:", N)
by = collections.Counter(d["exp"] for d in rows)
for k, v in by.most_common(): print(f"  growth {k:6s} {v:6d}  {100*v/N:5.1f}%")
print("\nby group:")
for g in sorted(set(d["group"] for d in rows)):
    c = collections.Counter(d["exp"] for d in rows if d["group"] == g); n = sum(c.values())
    print(" ", g, n, {k: v for k, v in c.most_common()})
print("\ntop law shapes per growth class:")
for k, _ in by.most_common():
    c = collections.Counter(d["form"] for d in rows if d["exp"] == k)
    print(f"  [{k}]", c.most_common(6))
cvs = np.array([d.get("cv", 0) for d in rows])
print("\ncv quantiles (non-bounded):", np.percentile(cvs[cvs > 0], [50, 90, 99, 100]))
