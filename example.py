"""Small example: Black's law of mixtures, T_f = (M1*T1 + M2*T2)/(M1 + M2), from a noisy
factorial experiment (3 levels per input, 81 rows, 1% noise on the measured T_f).

    python example.py
"""
from itertools import product

import numpy as np
import sympy as sp

from basin import discover

names = ["M1", "T1", "M2", "T2"]
levels = [[1.0, 2.0, 3.0], [10.0, 50.0, 90.0], [1.0, 2.0, 3.0], [20.0, 60.0, 100.0]]
X = np.array(list(product(*levels)))
M1, T1, M2, T2 = X.T
y_true = (M1 * T1 + M2 * T2) / (M1 + M2)
rng = np.random.default_rng(0)
y = y_true + rng.normal(0.0, 0.01 * np.abs(y_true))   # 1% noise, sd proportional to |y|

symbols = [sp.Symbol(v, positive=True) for v in names]
law, info = discover(X, y, symbols, details=True)
print("T_f =", sp.N(law, 4))
print("proposed powers:", info["proposed_powers"])
print(f"terms in the list: {info['n_terms_in_list']}, candidate laws scored: {info['n_candidates']}")
