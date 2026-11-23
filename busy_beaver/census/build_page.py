"""Build atlas.html from the classified censuses (BB(5) in census/, BB(3) and BB(4) in small/)."""
import json, collections, re, numpy as np
from fractions import Fraction

import os
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the busy_beaver/ folder
bb5 = json.load(open(f"{S}/census/classified.json"))
bb4 = json.load(open(f"{S}/small/bb4_classified.json"))
bb3 = json.load(open(f"{S}/small/bb3_classified.json"))
TEXT = json.load(open(f"{S}/census/class_text.json"))
META = json.load(open(f"{S}/census/page_meta.json"))
ORDER = ["0", "log", "1/3", "1/2", "1", "other"]
KEY = {"0": "bounded", "log": "log", "1/3": "third", "1/2": "half", "1": "one", "other": "other"}
CHAMP = dict(machine="1RB1LC_1RC1RB_1RD0LE_1LA1LD_---0LA", c=1.92019762050414, cv=1.48e-3)


def pretty(expr):
    def rnd(m):
        v = float(m.group(0))
        if abs(v) >= 1e4 or (abs(v) < 1e-3 and v != 0):
            return f"{v:.2e}"
        return f"{v:.3g}"
    e = re.sub(r"(?<![\w.])\d+\.\d+(e[-+]?\d+)?", rnd, expr)
    e = re.sub(r"\b1\.0?\*", "", e)
    e = e.replace("**", "^").replace("sqrt(t)", "√t").replace("*", "·").replace("u", "log₂t")
    return e


def form_pretty(f):
    return f.replace("**", "^").replace("sqrt(t)", "√t").replace("*", "·").replace("u", "log₂t")


def counts(rows):
    return dict(collections.Counter(KEY[d["exp"]] for d in rows))


groups = [
    dict(label="3 states", note="all non-halting machines", n=len(bb3), counts=counts(bb3)),
    dict(label="4 states", note="sample of 436,278", n=len(bb4), counts=counts(bb4)),
    dict(label="5 states", note="sample of 88,664,064", n=len(bb5), counts=counts(bb5)),
    dict(label="5 states, time limit", note="never used 12,289 cells", sub=True,
         n=sum(d["group"] == "time" for d in bb5), counts=counts([d for d in bb5 if d["group"] == "time"])),
    dict(label="5 states, space limit", note="used over 12,289 cells", sub=True,
         n=sum(d["group"] == "space" for d in bb5), counts=counts([d for d in bb5 if d["group"] == "space"])),
]


def hist_for(key, sub):
    r = np.array([d["rate"] for d in sub if d.get("rate") is not None and np.isfinite(d["rate"])])
    marker = None
    if len(r) < 20:
        return None
    if key == "bounded":
        r = np.clip(r, 1, 16); edges = np.arange(0.5, 16.6, 1)
        title = "Tape span in cells (bounded machines never leave it)"
        lab = lambda a, b: f"{int(a + 0.5)} cells"
        ticks = [(i, str(i + 1)) for i in range(0, 16, 3)]
    elif key == "log":
        edges = np.arange(0.5, 4.5001, 0.1); r = np.clip(r, 0.5, 4.4999)
        title = "Cells added each time t doubles"
        lab = lambda a, b: f"{a:.1f} to {b:.1f} cells per doubling"
        ticks = [(int(round((v - 0.5) / 0.1)), str(v)) for v in (1, 2, 3, 4)]
    elif key == "one":
        edges = np.arange(0, 1.0001, 0.01); r = np.clip(r, 0, 0.9999)
        title = "Speed of the growing edge, in cells per step"
        lab = lambda a, b: f"{a:.2f} to {b:.2f} cells/step"
        ticks = [(int(v * 100) - (1 if v == 1 else 0), s) for v, s in ((1/9, "1/9"), (1/5, "1/5"), (1/3, "1/3"), (1/2, "1/2"), (1, "1"))]
    elif key == "half":
        r = r ** 2; hi = 7.0
        edges = np.arange(0, hi + 1e-9, 0.05); r = np.clip(r, 0, hi - 1e-9)
        title = "c² for span ≈ c·√t; ticks mark whole numbers"
        lab = lambda a, b: f"c² from {a:.2f} to {b:.2f}"
        ticks = [(int(round(v / 0.05)), str(v)) for v in range(0, 7)]
        marker = dict(i=CHAMP["c"] ** 2 / 0.05, label="BB(5) champion, c² = 3.69")
    else:
        return None
    cnt, _ = np.histogram(r, edges)
    return dict(title=title, counts=cnt.tolist(), binLabels=[lab(a, b) for a, b in zip(edges, edges[1:])],
                ticks=ticks, marker=marker)


def pick_examples(key, sub, k=4):
    sub = sorted(sub, key=lambda d: d["i"])
    if key == "bounded":
        by = collections.defaultdict(list)
        for d in sub: by[d["span"]].append(d)
        spans = sorted(by)
        want = [4, 6, spans[len(spans) * 9 // 10], spans[-1]]
        out = [by[s][len(by[s]) // 2] for s in dict.fromkeys(w for w in want if w in by)]
    elif key in ("third", "other"):
        seen, out = set(), []
        for d in sub:
            if d["law"] not in seen: seen.add(d["law"]); out.append(d)
    else:
        by = collections.defaultdict(list)
        for d in sub:
            rr = d.get("rate")
            by[round(rr ** 2 if key == "half" else rr, 1 if key == "log" else 2)].append(d)
        common = sorted(by.items(), key=lambda kv: -len(kv[1]))
        out = [v[len(v) // 2] for _, v in common[:k]]
    res = []
    for d in out[:k]:
        law = f"span = {d['span']} cells" if key == "bounded" else "span ≈ " + pretty(d["law"])
        res.append(dict(machine=d["machine"], law=law))
    return res


def share(rows, e):
    return 100 * sum(d["exp"] == e for d in rows) / len(rows)


classes = []
for e in ORDER:
    key = KEY[e]
    sub = [d for d in bb5 if d["exp"] == e]
    T = TEXT[key]
    forms = collections.Counter(form_pretty(d["form"]) for d in sub).most_common(5) if key != "bounded" else []
    ex = pick_examples(key, sub) if sub else []
    if key == "third" and not sub:
        continue
    classes.append(dict(key=key, name=T["name"], color=T["color"], text=T["text"], steps=T["steps"],
                        across=[dict(n=n, v=share(rows, e)) for n, rows in ((3, bb3), (4, bb4), (5, bb5))],
                        forms=forms, hist=hist_for(key, sub), examples=ex))

# roughest sqrt fits next to the champion
half = sorted([d for d in bb5 if d["exp"] == "1/2"], key=lambda d: -d["cv"])
rough = [dict(machine=CHAMP["machine"], law=f"span ≈ 1.92·√t · fit error {CHAMP['cv']:.1e}", champ=True)]
rough += [dict(machine=d["machine"], law=f"span ≈ {pretty(d['law'])} · fit error {d['cv']:.1e}") for d in half[:7]]
cvs = np.array([d["cv"] for d in half])

data = dict(groups=groups, classes=classes, rough=dict(examples=rough, steps=TEXT["rough"]["steps"],
            text=TEXT["rough"]["text"].format(median=f"{np.median(cvs):.1e}", rank=int((cvs >= CHAMP["cv"]).sum()), n=len(cvs))),
            **META)
html = open(f"{S}/census/template.html").read().replace("/*DATA*/null", json.dumps(data, ensure_ascii=False))
open(f"{S}/census/atlas.html", "w").write(html)
print("wrote atlas.html", len(html) // 1024, "KB")
for g in groups: print(" ", g["label"], g["n"], g["counts"])
