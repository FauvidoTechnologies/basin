"""Figures for report/busy_beaver_report.tex. Run from busy_beaver/ with the report venv:

    report/.venv/bin/python report/figures.py            # all figures
    report/.venv/bin/python report/figures.py champion   # one figure

Reads results/ and pisot/ (pisot/long36.out for the convergence figure), runs short simulations
for the space-time diagrams, and calls bb6/simtxt (compile it first) for span curves.
"""
import gzip, math, os, subprocess, sys, collections, itertools
from fractions import Fraction
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__)); BB = os.path.dirname(HERE)
RES = os.path.join(BB, "results"); FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)
sys.path.insert(0, os.path.join(BB, "pisot"))

# reference palette (light mode, print): categorical slots in fixed order, ink, chrome
S1, S2, S3, S4, S5, S6, S7 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE500, SURF = "#256abf", "#fcfcfb"
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titleweight": "bold",
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
    "legend.frameon": False, "legend.fontsize": 8, "figure.dpi": 150, "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
W = 6.3  # text width in inches
TAPE = ListedColormap([SURF, BLUE500])


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf")); plt.close(fig); print("wrote", name)


def parse(m):
    T = {}
    for s, blk in enumerate(m.split("_")):
        for r in range(2):
            p = blk[3 * r:3 * r + 3]
            T[(s, r)] = None if p[2] not in "ABCDEFGH" else (int(p[0]), -1 if p[1] == "L" else 1, ord(p[2]) - 65)
    return T


def spacetime(m, steps, width=None, every=1, origin=None):
    """Rows = tape at every `every` steps. With width=None the window is every cell visited."""
    if width is None:
        T = parse(m); tape = collections.defaultdict(int); pos = lo = hi = 0; st = 0
        for t in range(steps):
            a = T[(st, tape[pos])]
            if a is None: break
            tape[pos] = a[0]; pos += a[1]; st = a[2]; lo = min(lo, pos); hi = max(hi, pos)
        origin, width = lo - 1, hi - lo + 3
    T = parse(m); tape = collections.defaultdict(int); pos = 0; st = 0; rows = []; heads = []
    o = -width // 2 if origin is None else origin
    for t in range(steps + 1):
        if t % every == 0:
            rows.append([tape[i] for i in range(o, o + width)]); heads.append(pos - o)
        a = T[(st, tape[pos])]
        if a is None: break
        tape[pos] = a[0]; pos += a[1]; st = a[2]
    return np.array(rows), np.array(heads)


def span_curve(machines):
    """span at t_k = round(2^(10+k/4)), k = 0..48, from bb6/simtxt."""
    lst = os.path.join(FIG, "_machines.txt"); out = os.path.join(FIG, "_span.bin")
    open(lst, "w").write("\n".join(machines) + "\n")
    subprocess.run([os.path.join(BB, "bb6", "simtxt"), lst, "22", out], check=True)
    K = 49; a = np.fromfile(out, dtype=np.int32).reshape(len(machines), 1 + 2 * K)
    os.remove(lst); os.remove(out)
    t = np.array([round(2 ** (10 + k / 4)) for k in range(K)], float)
    return t, a[:, 1:1 + K].astype(float)


def law(expr, t):
    u = np.log2(t)
    return eval(expr.replace("sqrt", "np.sqrt"), {"np": np, "t": t, "u": u}) * np.ones_like(t)


# --------------------------------------------------------------------------------------------
def classes():
    order = [("1", "linear"), ("0", "bounded"), ("1/2", r"$\sqrt{t}$"), ("log", r"$\log t$"),
             ("1/3", r"$t^{1/3}$"), ("other", "other")]
    sets = [("bb3_census.csv", "3 states (all 2,678)", S1, "o"), ("bb4_census.csv", "4 states (20,000)", S2, "s"),
            ("bb5_census.csv", "5 states (20,000)", S3, "D")]
    fig, ax = plt.subplots(figsize=(W, 2.5))
    for j, (f, lab, col, mk) in enumerate(sets):
        d = pd.read_csv(os.path.join(RES, f)); c = d.growth_class.value_counts(normalize=True) * 100
        xs, ys = [], []
        for i, (k, _) in enumerate(order):
            if c.get(k, 0) > 0: xs.append(c[k]); ys.append(i + (j - 1) * 0.18)
        ax.scatter(xs, ys, s=30, color=col, marker=mk, label=lab, zorder=3, edgecolor="white", linewidth=0.8)
    ax.set_xscale("log"); ax.set_xlim(0.003, 150)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([o[1] for o in order]); ax.invert_yaxis()
    ax.set_xlabel("share of machines (%, log scale)")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right")
    save(fig, "classes")


def span_examples():
    ex = [("linear, speed 1/3", "1RB0RC_0LC1RE_1LE1RD_1LD1RA_1RB---"),
          ("bounded (cycler)", "1RB---_0LC1RE_0LD1LD_0RE1LC_1RA1LB"),
          (r"$\sqrt{t}$ bouncer, $c=\sqrt{2}$", "1RB---_1LC0LB_0RD0LB_1RE1RC_1LD---"),
          (r"$\log t$ counter (binary)", "1RB1LA_1LC1RD_---0RB_0LA0RE_0LB1RD"),
          (r"$t^{1/3}$, $c^3=12$", "1RB1LA_1LA1RC_0RD0LB_0LA1RE_---1RD"),
          ("BB(5) champion (halts later)", "1RB1LC_1RC1RB_1RD0LE_1LA1LD_1RZ0LA")]
    cen = pd.read_csv(os.path.join(RES, "bb5_census.csv")).set_index("machine")
    b6 = pd.read_csv(os.path.join(RES, "bb6_fits.csv"))
    t, S = span_curve([m for _, m in ex])
    fig, axs = plt.subplots(2, 3, figsize=(W, 3.6))
    tt = np.logspace(12 * math.log10(2), 22 * math.log10(2), 200)
    for ax, (title, m), s in zip(axs.flat, ex, S):
        ax.axvspan(t[0], 2 ** 12, color=GRID, alpha=0.5, lw=0)
        ax.plot(t, s, "o", ms=3.4, mfc="white", mec=S1, mew=0.9, zorder=2, label="simulated span")
        expr = cen.law.get(m) if m in cen.index else b6[b6.machine == m].law.iloc[0]
        ax.plot(tt, law(str(expr), tt), color=S2, lw=1.3, zorder=3, label="basin law (fit window)")
        ax.set_xscale("log", base=2)
        if "counter" not in title: ax.set_yscale("log")
        ax.set_title(title, loc="left", fontsize=8.5)
        ax.set_xticks([2 ** 10, 2 ** 16, 2 ** 22]); ax.set_xticklabels(["$2^{10}$", "$2^{16}$", "$2^{22}$"])
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    for ax in axs[1]: ax.set_xlabel("steps $t$")
    for ax in axs[:, 0]: ax.set_ylabel("span (cells)")
    axs[0, 0].legend(loc="upper left", fontsize=7, handlelength=1.2)
    fig.tight_layout()
    save(fig, "span_examples")


def spacetime_fig():
    ex = [("translated cycler", "1RB0RC_0LC1RE_1LE1RD_1LD1RA_1RB---", 240, 1, None, None),
          ("cycler", "1RB---_0LC1RE_0LD1LD_0RE1LC_1RA1LB", 240, 1, None, None),
          ("bouncer", "1RB---_1LC0LB_0RD0LB_1RE1RC_1LD---", 12000, 12, None, None),
          ("BB(5) champion", "1RB1LC_1RC1RB_1RD0LE_1LA1LD_1RZ0LA", 60000, 40, None, None)]
    fig, axs = plt.subplots(1, 4, figsize=(W, 3.2))
    for ax, (title, m, n, every, w, o) in zip(axs, ex):
        R, H = spacetime(m, n, w, every, o); w = R.shape[1]
        ax.imshow(R, cmap=TAPE, aspect="auto", interpolation="nearest", extent=[0, w, n, 0])
        if every == 1: ax.plot(H + 0.5, np.arange(len(H)), "-", color=S2, lw=0.7)
        ax.set_title(title, loc="left", fontsize=8); ax.grid(False)
        ax.set_xticks([]); ax.tick_params(axis="y", labelsize=7)
        for sp in ax.spines.values(): sp.set_visible(True); sp.set_color(AXIS)
    axs[0].set_ylabel("step $t$ (down)")
    fig.tight_layout(w_pad=0.6)
    save(fig, "spacetime")


def counters_fig():
    """Tape at consecutive increments (head arriving at the counter's low end), low digit at left."""
    from weights2 import run as run2
    ex = [("binary (base 2)", "1RB---_0LC0LD_1RD0RB_1LE1RA_0RC0LE", ("L", 0)),
          (r"Lucas weights (base $\varphi$)", "1RB1LA_0LC0RD_---1RD_0LA1LE_0RC0RA", ("R", 0)),
          ("base 3/2, digits 0,1,2", "1RB1RA_1LC0RD_---1LD_0LE1LE_0RA0LD", ("L", 3))]
    fig, axs = plt.subplots(1, 3, figsize=(W, 3.0))
    for ax, (title, m, key) in zip(axs, ex):
        snaps = run2(m, 400_000)[key]
        snaps = snaps[len(snaps) // 2: len(snaps) // 2 + 120]
        w = max(len(c) for _, _, c in snaps)
        R = np.zeros((len(snaps), w))
        for i, (_, _, c) in enumerate(snaps):
            c = c if key[0] == "L" else c[::-1]; R[i, :len(c)] = c
        ax.imshow(R, cmap=TAPE, aspect="auto", interpolation="nearest")
        ax.set_title(title, loc="left", fontsize=8); ax.grid(False)
        ax.set_xlabel("cell (low digit at left)"); ax.tick_params(labelsize=7)
        for sp in ax.spines.values(): sp.set_visible(True); sp.set_color(AXIS)
    axs[0].set_ylabel("increment $n$ (120 consecutive)")
    fig.tight_layout(w_pad=0.8)
    save(fig, "counters_tape")


def speeds():
    d = pd.read_csv(os.path.join(RES, "bb5_census.csv")); v = d[d.growth_class == "1"].rate.astype(float)
    def frac(x):
        f = Fraction(x).limit_denominator(12); return f if abs(float(f) - x) < 2e-3 else None
    c = collections.Counter(frac(x) for x in v); c.pop(None, None)
    top = c.most_common(16)
    fig, ax = plt.subplots(figsize=(W, 2.4))
    for i, (f, n) in enumerate(top):
        even = (f.numerator + f.denominator) % 2 == 0
        ax.bar(i, n, width=0.72, color=S1 if even else S2, label=("period can be $b$ ($a+b$ even)" if even else "period at least $2b$ ($a+b$ odd)"))
    h, l = ax.get_legend_handles_labels(); u = dict(zip(l, h)); ax.legend(u.values(), u.keys(), loc="upper right")
    ax.set_xticks(range(len(top))); ax.set_xticklabels([f"{f.numerator}/{f.denominator}" if f.denominator > 1 else "1" for f, _ in top])
    ax.set_xlabel("speed $a/b$ of the growing edge (cells per step)"); ax.set_ylabel("machines")
    ax.grid(axis="x", visible=False)
    save(fig, "speeds")


def bouncers():
    d = pd.read_csv(os.path.join(RES, "bb5_census.csv")); c2 = d[d.growth_class == "1/2"].rate.astype(float) ** 2
    fig, ax = plt.subplots(figsize=(W, 2.3))
    bins = np.logspace(np.log10(0.3), np.log10(12), 260)
    ax.hist(c2, bins=bins, color=S1, lw=0)
    top = ax.get_ylim()[1] * 1.12; ax.set_ylim(0, top)
    for j, f in enumerate(["2/3", "1", "4/3", "3/2", "2", "8/3", "3", "4"]):
        x = float(Fraction(f)); ax.axvline(x, color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
        ax.text(x, top * (0.99 if j % 2 == 0 else 0.91), f, ha="center", va="top", fontsize=7, color=INK2,
                bbox=dict(fc="white", ec="none", pad=0.6))
    ax.set_xscale("log"); ax.set_xlim(0.3, 12)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(r"$c^2$ in span $\approx c\sqrt{t}$ (log scale)"); ax.set_ylabel("machines")
    save(fig, "bouncers")


def bases():
    d = pd.read_csv(os.path.join(RES, "bb5_counter_bases.csv.gz")); f = d.time_factor_per_new_cell
    lab = d.base.astype(str).copy()
    lab[np.abs(f / 1.5 ** .5 - 1) < 1e-3] = "sqrt(3/2)"
    lab[np.abs(f / 1.5 - 1) < 1e-3] = "3/2"
    lab[np.abs(f / (2 + 2 ** .5) ** .5 - 1) < 1e-3] = "sqrt(2+sqrt2)"
    lab[(f > 1.3603) & (f < 1.3613)] = "2^(4/9)?"
    names = {  # label: (display, kind)
        "sqrt2": (r"2, 2 cells/digit", "p"), "2": ("2, 1 cell/digit", "p"), "2^(1/3)": ("2, 3 cells/digit", "p"),
        "phi": (r"$\varphi$ (Fibonacci)", "p"), "2^(4/9)?": (r"$2^{4/9}$ (candidate, unchecked)", "c"),
        "2^(1/4)": ("2, 4 cells/digit", "p"), "sqrt3": ("3, 2 cells/digit", "p"), "2^(2/3)": ("2, 3 cells/2 digits", "p"),
        "3^(1/3)": ("3, 3 cells/digit", "p"), "sqrt(phi)": (r"$\varphi$, 2 cells/digit", "p"),
        "narayana": (r"Narayana $x^3=x^2+1$", "p"), "plastic": (r"plastic $x^3=x+1$", "p"),
        "sqrt(1+sqrt2)": (r"$1+\sqrt{2}$, 2 cells/digit", "p"), "sqrt(2+sqrt2)": (r"$2+\sqrt{2}$, 2 cells/digit", "p"),
        "plastic^2": ("plastic squared", "p"), "2^(1/6)": ("2, 6 cells/digit", "p"), "3^(1/4)": ("3, 4 cells/digit", "p"),
        "tribonacci": (r"tribonacci $x^3=x^2+x+1$", "p"), "sqrt(3/2)": ("3/2, 2 cells/digit", "r"),
        "3/2": (r"$\approx$1.5 per cell, irregular (1 behaviour)", "u"), "phi^(1/3)": (r"$\varphi$, 3 cells/digit", "p"),
        "2^(1/5)": ("2, 5 cells/digit", "p"), "unidentified": ("unidentified", "u")}
    c = lab.value_counts()
    rows = [(names[k][0], names[k][1], n) for k, n in c.items()]
    col = {"p": S1, "r": S2, "c": S7, "u": MUTED}
    fig, ax = plt.subplots(figsize=(W, 4.4))
    for i, (name, kind, n) in enumerate(rows):
        ax.barh(i, n, height=0.7, color=col[kind])
        ax.text(n * 1.08, i, f"{n:,}", va="center", fontsize=7, color=INK2)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=7.5); ax.invert_yaxis()
    ax.set_xscale("log"); ax.set_xlim(5, 1e6); ax.grid(axis="y", visible=False)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlabel("records among the 297,225 BB(5) counters (log scale)")
    from matplotlib.patches import Patch
    ax.legend([Patch(color=col[k]) for k in "prcu"],
              ["Pisot base", "rational base 3/2 (not Pisot)", "candidate", "unidentified"], loc="lower right")
    save(fig, "bases")


def weights_fig():
    from weights import run, fit
    from weights2 import run as run2, fit as fit2
    fig, ax = plt.subplots(figsize=(W, 2.6))
    L, R = run("1RB---_0LC0LD_1RD0RB_1LE1RA_0RC0LE", 1_500_000); w2, _, _ = fit(L, "L")
    L, R = run("1RB1LA_0LC0RD_---1RD_0LA1LE_0RC0RA", 1_500_000); wf, _, _ = fit(R, "R")
    S = run2("1RB1RA_1LC0RD_---1LD_0LE1LE_0RA0LD", 1_500_000); w3, _, _ = fit2(S[("L", 3)], "L")
    for w, lab, col, mk in [(w2, "binary machine", S1, "o"), (wf, r"$\varphi$ machine (Lucas numbers)", S3, "s"),
                            (w3, "base-3/2 machine", S2, "D")]:
        i = np.arange(len(w)); m = (np.abs(w) > 0.2) & (np.abs(w) < 1e3) & (i < 26)
        ax.plot(i[m], np.abs(w[m]), mk, ms=4.5, color=col, label=lab, mec="white", mew=0.6)
    ax.set_yscale("log"); ax.set_xlabel("cell index from the counter's low end"); ax.set_ylabel(r"fitted weight $|w_i|$")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.legend(loc="upper left")
    save(fig, "weights")


def champion():
    lines = open(os.path.join(RES, "champion_trace.txt")).read().splitlines()
    S = np.array([[int(x) for x in l.split()[1:]] for l in lines if l.startswith("S ")], float)
    t, r = S[:, 0], S[:, 1] / np.sqrt(S[:, 0])
    fig, ax = plt.subplots(figsize=(W, 2.5))
    ax.plot(t, r, color=S1, lw=1.2, label=r"simulated span$/\sqrt{t}$")
    for y, lab, ls in [(4 / 5 ** .5, r"$4/\sqrt{5}$ (at configurations)", (0, (3, 2))), (2, "2 (maximum in a step)", (0, (3, 2))),
                       (1.9236, "1.9236 (mean over $\\log t$)", "-"), (1.920, None, None)]:
        if lab is None: continue
        ax.axhline(y, color=MUTED if ls != "-" else S2, lw=0.9, ls=ls, zorder=1)
        ax.text(7e7, y, lab, fontsize=7, color=INK2, va="center", clip_on=False)
    ax.set_xscale("log", base=10); ax.set_xlim(1e3, 6e7); ax.set_ylim(1.72, 2.06)
    ax.axvline(47176870, color=INK2, lw=0.8); ax.text(4.3e7, 1.74, "halts at\n47,176,870", ha="right", fontsize=7, color=INK2)
    ax.set_xlabel("steps $t$ (log scale)"); ax.set_ylabel(r"span$/\sqrt{t}$")
    save(fig, "champion")


def roughness():
    b6 = pd.read_csv(os.path.join(RES, "bb6_fits.csv")); b5 = pd.read_csv(os.path.join(RES, "bb5_census.csv"))
    fig, ax = plt.subplots(figsize=(W, 2.6))
    sets = [(b5[b5.growth_class == "1/2"].cv_rel_mse, r"BB(5) non-halting, $\sqrt{t}$ class", S1),
            (b6[(b6.set == "random6") & (b6.growth_class == "1/2")].cv_rel_mse, r"random 6-state, $\sqrt{t}$ class", S3),
            (b6[(b6.set == "holdout") & (b6.growth_class == "1/2")].cv_rel_mse, r"BB(6) holdouts, $\sqrt{t}$ class", S2)]
    for x, lab, col in sets:
        x = np.sort(x[x > 0].astype(float)); ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=col, label=f"{lab} (n={len(x)})")
    cr = b6[(b6.set == "cryptid") & (b6.growth_class == "1/2")].cv_rel_mse.astype(float)
    ax.plot(cr, np.full(len(cr), -0.06), "|", color=INK, ms=9, mew=1.1, label=r"cryptids in the $\sqrt{t}$ class", clip_on=False)
    ax.set_xscale("log"); ax.set_ylim(-0.1, 1.02)
    ax.set_xlabel("roughness: basin's cross-validated relative error (log scale)"); ax.set_ylabel("cumulative share")
    ax.legend(loc="center", bbox_to_anchor=(0.47, 0.55), fontsize=7)
    save(fig, "roughness")


def convergence():
    path = os.path.join(BB, "pisot", "long36.out")
    if not os.path.exists(path): print("skip convergence (no pisot/long36.out)"); return
    phi = (1 + 5 ** .5) / 2
    groups = [("r32", "3/2, 2 cells/digit", 1.5 ** .5, S2, "D"), ("r3_2", r"$\approx$1.5 per cell, irregular", 1.5, S7, "v"),
              ("ctl_2", "2, 1 cell/digit", 2.0, S1, "o"), ("ctl_phi", r"$\varphi$", phi, S3, "s"),
              ("ctl_plastic", "plastic", 1.324717957244746, S4, "^")]
    Ls = [20, 22, 24, 26, 28, 30, 32, 34, 36]
    data = collections.defaultdict(list)
    for line in open(path):
        f = line.split(); g = f[0]; ts = np.array(f[3:], dtype=float); data[g].append(ts)
    fig, ax = plt.subplots(figsize=(W, 2.7))
    for g, lab, target, col, mk in groups:
        if g not in data: continue
        med = []
        for L in Ls:
            errs = []
            for ts in data[g]:
                t = ts[ts <= 2 ** L]; late = t[t > 2 ** 14]; q = (len(late) - 1) // 6 * 6
                if q >= 6: errs.append(abs((t[-1] / t[-1 - q]) ** (1 / q) / target - 1))
            med.append(np.median(errs) if errs else np.nan)
        ax.plot(Ls, med, "-" + mk, color=col, ms=4.5, mec="white", mew=0.6, label=f"{lab} (n={len(data[g])})")
    ax.set_yscale("log"); ax.set_xticks(Ls); ax.set_xticklabels([f"$2^{{{L}}}$" for L in Ls])
    ax.set_xlabel("run length (steps)"); ax.set_ylabel("median relative error\nof the time factor per cell")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=7, ncol=3)
    save(fig, "convergence")


def wobble():
    """Detrended ln(span) against ln t for four machines, and the multipliers of all detections."""
    W_ = os.path.join(BB, "wobble")
    spans = {l.split()[0]: np.array(l.split()[2:], float) for l in open(os.path.join(W_, "dense32.out"))}
    ex = [("Antihydra: $\\lambda = 3/2$", "1RB1RA_0LC1LE_1LD1LC_1LA0LB_1LF1RE_---0RA", S1),
          ("cryptid 6: $\\lambda \\approx 4/3$", "1RB1LD_1RC1RE_0LA1LB_0LD1LC_1RF0RA_---0RC", S3),
          ("cryptid 20: $\\lambda = 8/3$", None, S7),
          ("ordinary BB(5) bouncer", None, MUTED)]
    sp = pd.read_csv(os.path.join(BB, "bb6", "special.csv")).set_index("name").machine
    ex[1] = (ex[1][0], sp["cryptid 6"], S3); ex[2] = (ex[2][0], sp["cryptid 20"], S7)
    wr = pd.read_csv(os.path.join(W_, "wobble_results.csv"))
    ex[3] = (ex[3][0], wr[wr.set == "bb5_sqrt"].sort_values("R2").iloc[100].machine, MUTED)
    fig = plt.figure(figsize=(W, 4.6)); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.15])
    for k, (title, m, col) in enumerate(ex):
        ax = fig.add_subplot(gs[k // 2, k % 2]); s = spans[m]
        j = np.arange(len(s)); t = 2.0 ** (10 + j / 64); ok = (s > 0) & (t >= 2 ** 12)
        x, y = np.log(t[ok]), np.log(s[ok]); pp = np.polyfit(x, y, 1); r = y - np.polyval(pp, x)
        ax.plot(t[ok], r, color=col, lw=1.0)
        ax.set_xscale("log", base=2); ax.set_xticks([2 ** 12, 2 ** 22, 2 ** 32]); ax.set_xticklabels(["$2^{12}$", "$2^{22}$", "$2^{32}$"])
        ax.set_title(title, loc="left", fontsize=8); ax.set_ylim(-0.16, 0.12); ax.tick_params(labelsize=7)
        if k % 2 == 0: ax.set_ylabel("residual of\nln span", fontsize=7.5)
    ax = fig.add_subplot(gs[2, :])
    rf = pd.read_csv(os.path.join(W_, "refined.csv")).drop_duplicates("machine")
    groups = [("holdout", "BB(6) holdouts", S2, "o"), ("cryptid", "named cryptids", S1, "D"),
              ("halter5", "BB(5) champion family", S3, "s"), ("bb5_sqrt", "BB(5) non-halting", INK2, "^")]
    rng = np.random.default_rng(0)
    for i, (g, lab, col, mk) in enumerate(groups):
        v = rf[rf.set == g].lam_minima.dropna(); v = v[v < 4.5]
        ax.scatter(v, i + rng.uniform(-0.18, 0.18, len(v)), s=16, color=col, marker=mk, edgecolor="white", linewidth=0.5, zorder=3)
    for f in ["4/3", "3/2", "5/3", "2", "8/3", "10/3", "4"]:
        x = float(Fraction(f)); ax.axvline(x, color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
        ax.text(x, 3.65, f, ha="center", va="bottom", fontsize=7, color=INK2)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels([g[1] for g in groups], fontsize=7.5)
    ax.set_ylim(-0.6, 3.6); ax.grid(axis="y", visible=False); ax.set_xscale("log")
    ax.set_xticks([1.25, 1.5, 2, 2.5, 3, 4]); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"multiplier $\lambda$: span ratio between successive wobble minima (log scale)")
    fig.tight_layout(h_pad=0.8)
    save(fig, "wobble")


def rulecheck():
    """basin's roughness against the scale-invariance index for the BB(5) census sqrt(t) machines,
    coloured by the rule read from the tape."""
    d = pd.read_csv(os.path.join(RES, "bb5_census_sqrt_scaleinv.csv"))
    fig, a = plt.subplots(figsize=(W * 0.72, 2.8))
    fades = d[d.cls == "fades"]
    a.scatter(fades.S, fades.cv_rel_mse, s=8, color=MUTED, alpha=0.55, linewidth=0, label=f"wobble fades ({len(fades)})")
    for key, lab, col, mk in [("affine x2", "doubler, one exact formula", S1, "o"),
                              ("residue", "Collatz-like, exact formula per residue", S2, "D")]:
        v = d[d.tape_rule.fillna("").str.startswith(key)]
        a.scatter(v.S, v.cv_rel_mse, s=32, color=col, marker=mk, edgecolor="white", linewidth=0.6, zorder=3,
                  label=f"{lab} ({len(v)})")
    a.set_yscale("log"); a.set_xlabel("scale-invariance index $S$")
    a.set_ylabel("basin roughness (cv rel. MSE)")
    a.axvline(0.25, color=AXIS, lw=0.8, ls=(0, (2, 2))); a.set_xlim(-0.15, 1.35)
    a.legend(loc="upper left", bbox_to_anchor=(0.3, 1.0), fontsize=7, handletextpad=0.3, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, "rulecheck")


def screenrules():
    """Distinct scale-invariant BB(5) behaviours by the rule read from the tape."""
    r = pd.read_csv(os.path.join(RES, "bb5_screen_tape_rules.csv"))
    u = pd.read_csv(os.path.join(RES, "bb5_screen_full_scaleinv.csv")).drop_duplicates("fp")
    lam = r.frac.astype(str).replace({"9/4": "3/2"})
    rows = []
    for f in ["3/2", "4/3", "5/2", "5/4", "5/3", "7/4", "8/3"]:
        rows.append((f"$\\times{f}$", ((r.cat == "Collatz-like (fractional)") & (lam == f)).sum(), S2))
    for f in ["2", "3", "4"]:
        rows.append((f"$\\times{f}$", ((r.cat == "counter (one formula)") & (lam == f)).sum(), S1))
    rows.append(("$\\times5$, $\\times6$", ((r.cat == "counter (one formula)") & lam.isin(["5", "6"])).sum(), S1))
    rows.append(("teeth 2", (u.frac.astype(str) == "2").sum(), MUTED))
    rows.append(("no fraction", (u.frac.astype(str) == "none").sum(), MUTED))
    fig, ax = plt.subplots(figsize=(W * 0.8, 3.3))
    y = np.arange(len(rows))[::-1]
    ax.barh(y, [v for _, v, _ in rows], color=[c for _, _, c in rows], height=0.62, edgecolor=SURF, linewidth=1)
    for yi, (lab, v, c) in zip(y, rows):
        ax.text(v * 1.12, yi, f"{v:,}", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([l for l, _, _ in rows], fontsize=7.5)
    ax.set_xscale("log"); ax.set_xlim(0.8, 1.5e6); ax.grid(axis="y", visible=False)
    ax.set_xlabel("distinct behaviours (log scale)")
    for yi, lab, c in [(y[3], "Collatz-like,\nby residue", S2), (y[9], "exponential\ncounters", S1), ((y[11] + y[12]) / 2, "not read\non the tape", MUTED)]:
        ax.text(6.0e4, yi, lab, va="center", ha="left", fontsize=7.5, color=INK2, fontweight="bold")
        ax.plot([4.2e4], [yi], marker="s", ms=7, color=c, clip_on=False)
    fig.tight_layout()
    save(fig, "screenrules")


FIGS = dict(classes=classes, span_examples=span_examples, spacetime=spacetime_fig, counters=counters_fig,
            speeds=speeds, bouncers=bouncers, bases=bases, weights=weights_fig, champion=champion,
            roughness=roughness, convergence=convergence, wobble=wobble,
            rulecheck=rulecheck, screenrules=screenrules)
if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGS): FIGS[name]()
