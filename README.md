# basin

**Find the formula behind a table of measurements.** Give basin some inputs
`x₁ … xₙ` and one measured output `y`, and it returns a readable law
`y = N(x)/D(x)` as a sympy expression. For example, from noisy mixing data it recovers
Black's law `(M1·T1 + M2·T2)/(M1 + M2)`.

It belongs to the same family as [PySR](https://github.com/MilesCranmer/PySR) and
[PySINDy](https://github.com/dynamicslab/pysindy), and borrows from both sides:

- **like BACON** (Langley's discovery system), it proposes terms from the data itself,
  including fractional powers such as `x^(3/2)`;
- **like SINDy**, it fits sparse sums of those terms, as numerator and denominator;
- it **picks between candidate laws by cross-validation**, preferring the simplest law
  whose held-out error is within one standard error of the best.

On the AI Feynman benchmark it is level with PySR on algebraic laws, at about a fifth of
the CPU time. It cannot write `sin`, `cos`, `exp` or `log` ([benchmarks](#benchmarks)).

## Install

```bash
pip install basin-sr          # imported as `basin`; Python 3.10+
```

## Quick start

```python
from itertools import product

import numpy as np
import basin

# Black's law of mixtures from a 3-level factorial experiment (81 rows) with 1% noise
X = np.array(list(product([1, 2, 3], [10, 50, 90], [1, 2, 3], [20, 60, 100])), float)
M1, T1, M2, T2 = X.T
y = (M1 * T1 + M2 * T2) / (M1 + M2)
y = y + np.random.default_rng(0).normal(0, 0.01 * np.abs(y))

law = basin.discover(X, y, ["M1", "T1", "M2", "T2"])
print(law.evalf(4))    # (1.004*M1*T1 + 1.002*M2*T2)/(1.0*M1 + 1.003*M2), about 10 s
```

`example.py` in this repo runs the same thing.

## API

```python
law = basin.discover(X, y, symbols=None, details=False)
```

- `X`: array of shape `(rows, inputs)`; `y`: the measured output, shape `(rows,)`.
- `symbols`: input names (strings or sympy symbols; strings become positive symbols).
  The default is `x0, x1, …`.
- Returns a sympy expression, or `None` if no candidate law could be fitted.
- `details=True` also returns a dict: the proposed powers and their estimates, the
  number of terms and candidate laws, and the chosen law's cross-validated error.

```python
with basin.settings(POWER_FREE=True):       # change switches for one block only
    law = basin.discover(X, y, names)
```

`basin.settings` takes any of the settings in the table below and restores the old
values on exit. The switches are module globals, so this is not thread-safe.
`basin.propose_powers` and `basin.term_list` expose steps 1 and 2 on their own.

## How it works

1. **Propose powers.** For each input, basin takes the rows where every other input is
   held fixed (BACON's "slices" of a factorial experiment) and fits a power law
   `log|y| ~ p·log x` across all slices at once. Simple fractions near `p` become
   candidate terms (denominator ≤ 3, |p| ≤ 3, within 3 standard errors), and so does
   their product. Without slices, as in a time series, it uses one global log fit.
2. **Build a term list:** `1`, each input, each product of two inputs, plus the
   proposed terms.
3. **Fit fractions as a linear problem** (as in SINDy-PI): `y = N/D` becomes
   `y·D = N`, with `N` and `D` sparse sums of terms. Candidates are every law with up
   to 4 terms, plus larger ones (up to 5 coefficients) from STLSQ and forward
   selection.
4. **Decide by held-out error:** 5-fold cross-validation repeated 5 times on relative
   error. Pick the fewest terms within one standard error of the best; ties go to the
   lower total power.

**What it cannot write:** `sin`, `cos`, `exp`, `log`; powers with a denominator above 3
or above 3 in size; laws with more than 5 coefficients; products of more than two inputs
unless step 1 proposes them. It also cannot write fractional powers of inputs that are
not all positive. Step 1 is skipped if any `y` is exactly 0.

| setting | default | meaning |
|---|---|---|
| `MAX_FREE` | 5 | most free coefficients in a law (N and D together) |
| `EXHAUSTIVE_TERMS` | 4 | every law with up to this many terms is tried |
| `SCREEN_KEEP` | 30 | per law size, candidates kept after a quick screen |
| `POWER_SE`, `POWER_FLOOR` | 3, 0.03 | step 1 tolerance: 3 standard errors, at least 0.03 |
| `POWER_KEEP` | 2 | at most 2 proposed fractions per input |
| `FOLDS`, `REPEATS`, `CV_SEED` | 5, 5, 0 | cross-validation (leave-one-out below 10 rows) |
| `CV_FLOOR` | 1e-10 | errors closer than this count as equal |
| `POWER_DENOMINATORS` | (1, 2, 3) | step 1's fractions k/d use these d |
| `POWER_FREE`, `POWER_FREE_DIGITS` | False, 2 | if no fraction fits, propose the estimate rounded to 2 decimals |

`POWER_DENOMINATORS` and `POWER_FREE` were added after the Feynman benchmarks below. At
their defaults, basin gives exactly the benchmarked results.

## Benchmarks

basin was developed under the name "the hybrid", and the result files, the CSVs and
`findings/` still call it `Hybrid`. For each test, the rules were written down before it
ran, and the method was not changed after seeing any result. The exceptions are marked:
Tests 13–14 (it was designed after seeing these) and the real-data runs. "Exact law"
means the right terms with every number within 5% of its true value (equivalent
rearrangements count). Noise is Gaussian on the measured output, with sd = level × |y|.

### Summary: basin against PySR

- **Algebraic laws: level.** On the 59 laws without `sin`, `cos`, `exp` or `arcsin`, the
  two find about the same laws at every noise level (after-the-run split, not
  significant).
  - Basin is more reliable on some laws: 10 of 10 runs at every noise level, where
    PySR gets 3–4 of 10 on I.18.4 `(m1*r1 + m2*r2)/(m1 + m2)` and 2–5 on III.21.20
    `-rho_c_0*q*A_vec/m`.
  - PySR finds II.37.1 `mom*(1 + chi)*B` in 10 of 10 runs at every noise level, and
    I.16.6 on clean data. Basin finds neither in any run. `chi` gets no proposed
    power, so the term `mom*chi*B` is not in its list.
- **Trig and exponential laws: PySR only.** Basin cannot write them. This is the
  whole of PySR's lead in the overall counts.
- **Noise: basin's law count barely moves** (43, 43, 42, 41 laws). But its
  coefficients shrink with noise (Test 15), so PySR's answers are closer to the true law
  at 5–10% noise.
- **Clean data: PySR sometimes adds a tiny extra term**, e.g.
  `omega_0*v/(c - v + 7.18e-8) + omega_0`. That costs it 11 laws that it finds at 1%
  noise.
- **Small designs (5–27 rows, Tests 13–14): basin holds up better.** PySR builds
  long formulas that fit the few points. Neither finds the ideal gas law, because
  `T + 273` and its first-order twin fit equally well over the narrow range of `T`.
- **Fewer rows on the 79 Feynman laws (Tests 16–17): the result depends on how PySR
  chooses its formula.**
  - With PySR's default (`"best"`, Test 16) basin beats it at every noise level, at
    about 100 and about 30 rows (every p ≤ 0.003). PySR finds 4, 2 and 1 laws with noise
    at about 100 rows, and none at about 30.
  - With `"score"` (Test 17) the two are level at about 100 rows. At about 30 rows PySR
    beats basin with noise (41 laws to 10, 41 to 10, 36 to 9), mostly on the
    4-input laws, which then have 2 levels per input.
  - So Test 16 does not confirm Tests 13–14's "basin holds up better with few
    rows" against PySR at its best. Tests 13–14 used PySR's default choice too.
- **Cost: basin needs about a fifth of PySR's CPU time** (34 against 181 CPU-hours
  for Test 15, 13 against 55 at about 100 rows), and about a quarter to a third at about
  30 rows (12 against 43).

### Test 15: the AI Feynman benchmark (79 laws, basin was not designed on these)

- **Laws:** every AI Feynman equation (Udrescu & Tegmark 2020) with at most 4 inputs, 79
  laws. No other filtering.
- **Data:** a full factorial grid over the benchmark's ranges, about 1000 rows (1 input:
  1000 levels; 2 inputs: 32 per input; 3: 10; 4: 6).
- **Runs:** no noise (1 run per law), 1%, 5% and 10% noise (10 runs per law each).
- **PySR 2.7.2:** defaults plus every function the 79 laws use (`+ - * /`, `sqrt`, `exp`,
  `sin`, `cos`, `arcsin`), maximum size 40 (the largest true formula needs about 34),
  100 iterations, deterministic. PySR's own "best" model selection.
- **Scoring:** the answer's numbers are re-fitted to the true law at 500 random points.
  Right terms means the re-fitted formula reproduces the law exactly and no number goes
  to zero. Exact law means right terms and every number within 5% of its re-fitted value.
  This uses the fixed coefficient check from `rescore.py`.

**Laws found in more than half the runs, all five methods:**

| method | no noise (1 run per law) | 1% noise (10 runs) | 5% noise (10 runs) | 10% noise (10 runs) |
|---|---|---|---|---|
| BACON (original) | 40 of 79 | 28 of 79 | 19 of 79 | 14 of 79 |
| BACON-robust | 44 of 79 | 41 of 79 | 36 of 79 | 14 of 79 |
| PySINDy | 35 of 79 | 0 of 79 | 0 of 79 | 0 of 79 |
| **basin** | 43 of 79 | 43 of 79 | 42 of 79 | 41 of 79 |
| **PySR** | 42 of 79 | 50 of 79 | 50 of 79 | 49 of 79 |

**The same, split by kind of law.** The split was made after the run, from the true
formulas:
- **A:** one product of powers, e.g. `q1*q2/(4*pi*epsilon*r**2)`.
- **B:** other algebraic laws: sums, fractions with sums, roots of sums.
- **C:** laws with `sin`, `cos`, `exp` or `arcsin`.

| kind of law | noise | BACON (original) | BACON-robust | PySINDy | basin | PySR |
|---|---|---|---|---|---|---|
| A: one product of powers (35 laws) | no noise | 34 | 35 | 32 | 35 | 30 |
|  | 1% noise | 28 | 34 | 0 | 35 | 34 |
|  | 5% noise | 18 | 34 | 0 | 35 | 34 |
|  | 10% noise | 13 | 12 | 0 | 35 | 34 |
| B: other algebraic (24 laws) | no noise | 6 | 9 | 3 | 8 | 6 |
|  | 1% noise | 0 | 7 | 0 | 8 | 7 |
|  | 5% noise | 1 | 2 | 0 | 7 | 7 |
|  | 10% noise | 1 | 2 | 0 | 6 | 7 |
| C: sin, cos, exp, arcsin (20 laws) | no noise | 0 | 0 | 0 | 0 | 6 |
|  | 1% noise | 0 | 0 | 0 | 0 | 9 |
|  | 5% noise | 0 | 0 | 0 | 0 | 9 |
|  | 10% noise | 0 | 0 | 0 | 0 | 8 |

**basin against PySR, law by law.** Each count is the number of laws where that method
finds the exact law in more runs than the other. The sign test is two-sided. The "all 79"
rows were the test's decision rule. The splits by kind were made after the run, so they
describe the results rather than test them.

| laws | noise | basin finds the law in more runs | PySR finds it in more runs | same | sign test |
|---|---|---|---|---|---|
| all 79 | no noise | 8 | 7 | 64 | p = 1 |
|  | 1% noise | 6 | 21 | 52 | p = 0.0059 |
|  | 5% noise | 8 | 21 | 50 | p = 0.024 |
|  | 10% noise | 6 | 15 | 58 | p = 0.078 |
| A: one product of powers (35) | no noise | 5 | 0 | 30 | p = 0.062 |
|  | 1% noise | 3 | 0 | 32 | p = 0.25 |
|  | 5% noise | 5 | 0 | 30 | p = 0.062 |
|  | 10% noise | 3 | 0 | 32 | p = 0.25 |
| B: other algebraic (24) | no noise | 3 | 1 | 20 | p = 0.62 |
|  | 1% noise | 3 | 5 | 16 | p = 0.73 |
|  | 5% noise | 3 | 4 | 17 | p = 1 |
|  | 10% noise | 3 | 3 | 18 | p = 1 |
| C: sin, cos, exp, arcsin (20) | no noise | 0 | 6 | 14 | p = 0.031 |
|  | 1% noise | 0 | 16 | 4 | p = 3.1e-05 |
|  | 5% noise | 0 | 17 | 3 | p = 1.5e-05 |
|  | 10% noise | 0 | 12 | 8 | p = 0.00049 |
| A + B: the 59 algebraic laws | no noise | 8 | 1 | 50 | p = 0.039 |
|  | 1% noise | 6 | 5 | 48 | p = 1 |
|  | 5% noise | 8 | 4 | 47 | p = 0.39 |
|  | 10% noise | 6 | 3 | 50 | p = 0.51 |

**How close the answers predict, against the true law.** These counts include runs that
did not find the exact law. Error is the relative RMS error at 200 random points in the
ranges. Each cell counts the runs whose answer is within that error.

| laws | noise | runs | within 0.1%: basin | PySR | within 1%: basin | PySR | within 5%: basin | PySR |
|---|---|---|---|---|---|---|---|---|
| all 79 | no noise | 79 each | 43 | 66 | 53 | 73 | 61 | 75 |
|  | 1% noise | 790 each | 423 | 483 | 528 | 718 | 618 | 766 |
|  | 5% noise | 790 each | 0 | 261 | 483 | 554 | 580 | 749 |
|  | 10% noise | 790 each | 0 | 214 | 0 | 470 | 592 | 675 |
| 59 algebraic | no noise | 59 each | 43 | 53 | 50 | 57 | 55 | 57 |
|  | 1% noise | 590 each | 423 | 377 | 500 | 557 | 558 | 580 |
|  | 5% noise | 590 each | 0 | 206 | 472 | 429 | 540 | 570 |
|  | 10% noise | 590 each | 0 | 176 | 0 | 380 | 543 | 523 |

**Basin's numbers shrink with noise.** On runs that find the exact law, the median
prediction error is:

| noise | basin: exact runs | median error | PySR: exact runs | median error |
|---|---|---|---|---|
| no noise | 43 | 0.00% | 42 | 0.00% |
| 1% noise | 430 | 0.04% | 522 | 0.02% |
| 5% noise | 422 | 0.57% | 517 | 0.09% |
| 10% noise | 410 | 2.27% | 491 | 0.21% |

On simple products basin answers `0.9927*Nn*mu` at 5% noise and `0.9751*Nn*mu` at
10% (`F = mu*Nn`; the same on `m*g*z`, `q/C`, `Ef*q2`), in every run. PySR answers `Nn*mu`.

- **The cause is basin's fitting rule.** Step 4 fits on error relative to the noisy
  `y`, `(prediction − y)/y`. Dividing by a noisy `y` pulls every coefficient down by about
  2σ²: 0.5% at 5% noise, 2% at 10%.
- **The verdicts don't change.** These answers are still within the 5% check, so this
  does not change any count above. But at around 15–20% noise the bias alone would push
  coefficients past 5%.
- **Part of PySR's edge here is the benchmark.** Many Feynman laws have the coefficient
  1, which PySR can write with no fitted number at all.

**Run time** (all methods running 18 jobs in parallel on a 22-thread laptop at about
1.2 GHz. Rerun alone on this laptop, four basin runs took 3 to 9 times less time):

| inputs | laws | basin, median seconds per run | PySR, median seconds per run |
|---|---|---|---|
| 1 | 1 | 6 | 131 |
| 2 | 15 | 17 | 102 |
| 3 | 36 | 43 | 362 |
| 4 | 27 | 69 | 149 |
| all | 79 | 46 (total 34 CPU-hours) | 218 (total 181 CPU-hours) |

### Test 16: Test 15 with fewer rows (about 100 and about 30)

- **Data:** the same 79 laws and the same kind of grid, with fewer levels per input.
  - About 100 rows: 1 input 100 levels, 2 inputs 10, 3 inputs 5, 4 inputs 3 (81 rows).
  - About 30 rows: 30, 5, 3 and 2 levels (16 rows for 4 inputs).
- **Everything else as in Test 15:** noise levels, seeds, settings and scoring. The rules
  were written down before running, and no method was changed.

**Laws (of 79) whose exact law the method finds in more than half its runs** (fixed
coefficient check):

| method | rows per law | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|---|
| BACON (original) | about 1000 (Test 15) | 40 | 28 | 19 | 14 |
|  | about 100 | 40 | 28 | 17 | 17 |
|  | about 30 | 29 | 20 | 15 | 10 |
| BACON-robust | about 1000 (Test 15) | 44 | 41 | 36 | 14 |
|  | about 100 | 33 | 27 | 14 | 0 |
|  | about 30 | 14 | 11 | 0 | 0 |
| PySINDy | about 1000 (Test 15) | 35 | 0 | 0 | 0 |
|  | about 100 | 34 | 1 | 0 | 0 |
|  | about 30 | 30 | 6 | 0 | 0 |
| basin | about 1000 (Test 15) | 43 | 43 | 42 | 41 |
|  | about 100 | 43 | 43 | 41 | 37 |
|  | about 30 | 43 | 29 | 26 | 23 |
| PySR | about 1000 (Test 15) | 42 | 50 | 50 | 49 |
|  | about 100 | 24 | 4 | 2 | 1 |
|  | about 30 | 28 | 0 | 0 | 0 |

**basin against PySR, head to head.** Each cell: laws where basin finds the exact
law in more runs than PySR, then laws where PySR does, then the two-sided sign test.

| rows per law | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|
| about 1000 (Test 15) | 8 to 7, p = 1 | 6 to 21, p = 0.0059 | 8 to 21, p = 0.024 | 6 to 15, p = 0.078 |
| about 100 | 20 to 1, p = 2.1e-05 | 43 to 15, p = 0.00031 | 42 to 11, p = 2.2e-05 | 38 to 10, p = 6.2e-05 |
| about 30 | 19 to 4, p = 0.0026 | 43 to 2, p = 5.9e-11 | 41 to 1, p = 2e-11 | 36 to 0, p = 2.9e-11 |

With the original coefficient check, the 1% verdict at about 100 rows is level (36 to 21,
p = 0.063). Every other cell above is significant either way.

**Basin by number of inputs** (laws found in more than half the runs; no noise, 1%,
5%, 10%):

| inputs (laws) | Test 15 | about 100 rows | about 30 rows |
|---|---|---|---|
| 1 (1) | 1000 levels: 0, 0, 0, 0 | 100 levels: 0, 0, 0, 0 | 30 levels: 0, 0, 0, 0 |
| 2 (15) | 32 levels: 13, 13, 12, 12 | 10 levels: 13, 13, 12, 11 | 5 levels: 13, 12, 11, 11 |
| 3 (36) | 10 levels: 17, 17, 17, 16 | 5 levels: 17, 17, 17, 15 | 3 levels: 17, 14, 13, 12 |
| 4 (27) | 6 levels: 13, 13, 13, 13 | 3 levels: 13, 13, 12, 11 | 2 levels: 13, 3, 2, 0 |

**Median prediction error** (normalized RMS error at 200 random points, all runs; both
methods return a formula in every run). basin / PySR:

| rows per law | 1% noise | 5% noise | 10% noise |
|---|---|---|---|
| about 1000 (Test 15) | 0.09% / 0.06% | 0.77% / 0.29% | 2.6% / 0.60% |
| about 100 | 0.20% / 0.70% | 1.3% / 3.4% | 3.4% / 6.7% |
| about 30 | 0.76% / 2.0% | 3.4% / 8.4% | 6.7% / 15% |

- **PySR, with its default choice, collapses with noise and few rows; basin does
  not** (with `"score"` it does not collapse either: Test 17). At about 100 rows,
  basin is almost unchanged from Test 15. PySR finds 4, 2 and 1 laws with noise, and none
  at about 30 rows.
- **Why PySR collapses** (checked on one run, I.12.1 `F = mu*Nn`, about 100 rows, 1%
  noise):
  - PySR keeps `Nn*mu` among its formulas, with the highest score of any.
  - Its default choice (`model_selection="best"`) only considers formulas whose loss is
    within 1.5 times the lowest. With 100 rows, long formulas full of `cos` and `exp` fit
    the noise down to a loss of 0.0057, while `Nn*mu` has 0.0135, the noise level.
  - So it answers `Nn*(-0.002042*Nn*mu*(cos(Nn*exp(mu)) - cos(Nn*mu*cos(mu))) + mu)`.
  - With Test 15's 1024 rows, long formulas reach only 0.0093 against 0.0098, so
    `Nn*mu` is chosen.
  - Test 17 (below) reran these PySR runs with its `"score"` choice.
- **Basin's weak spot is 2 levels per input with noise:** the 4-input laws at about
  30 rows (16 rows) drop from 13 to 3, 2 and 0. Two points cannot tell a straight line
  from a curve, so many wrong formulas fit the data equally well. It also loses the other
  algebraic laws (kind B) with noise at about 30 rows: 6, 0 and 0, against 8, 7 and 6 in
  Test 15.

### Test 17: Test 16 again, PySR choosing with `model_selection="score"`

- **Setup:** Test 16's data, seeds and PySR settings, with one change. PySR picks the
  formula with the highest score among the ones it kept (`"score"`), with no loss
  cut-off. Its default `"best"` only looks at formulas within 1.5 times the lowest loss.
  Basin was not rerun: its Test 16 answers are compared run for run.
- **Check:** PySR's search repeated exactly. In all 4,898 runs the formula `"best"`
  would pick is Test 16's answer, so only the choice changed.
- Rules were written down before running. PySR got
  this second setting after its first failed; basin got no second setting.

**Laws found in more than half the runs** (of 79; fixed coefficient check):

| method | rows per law | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|---|
| basin | about 100 | 43 | 43 | 41 | 37 |
| PySR, `"best"` (Test 16) | about 100 | 24 | 4 | 2 | 1 |
| PySR, `"score"` | about 100 | 50 | 45 | 38 | 34 |
| basin | about 30 | 43 | 29 | 26 | 23 |
| PySR, `"best"` (Test 16) | about 30 | 28 | 0 | 0 | 0 |
| PySR, `"score"` | about 30 | 45 | 46 | 37 | 29 |

**basin against PySR with `"score"`, head to head.** Each cell: laws where basin
finds the exact law in more runs, then laws where PySR does, then the two-sided sign test.

| rows per law | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|
| about 100 | 3 to 10, p = 0.092 | 10 to 14, p = 0.54 | 14 to 15, p = 1 | 14 to 15, p = 1 |
| about 30 | 8 to 10, p = 0.81 | 10 to 41, p = 1.5e-05 | 10 to 41, p = 1.5e-05 | 9 to 36, p = 6.6e-05 |

**By kind of law** (A one product of powers, B other algebraic, C with `sin`, `cos`, `exp`
or `arcsin`), laws found in more than half the runs, basin / PySR with `"score"`:

| kind | rows per law | no noise | 1% noise | 5% noise | 10% noise |
|---|---|---|---|---|---|
| A (35) | about 100 | 35 / 34 | 35 / 32 | 35 / 28 | 35 / 25 |
| A (35) | about 30 | 35 / 32 | 23 / 31 | 26 / 26 | 23 / 23 |
| B (24) | about 100 | 8 / 8 | 8 / 7 | 6 / 6 | 2 / 5 |
| B (24) | about 30 | 8 / 5 | 6 / 7 | 0 / 5 | 0 / 2 |
| C (20) | about 100 | 0 / 8 | 0 / 6 | 0 / 4 | 0 / 4 |
| C (20) | about 30 | 0 / 8 | 0 / 8 | 0 / 6 | 0 / 4 |

**Median prediction error** at 200 random points, basin / PySR with `"score"` (PySR with
`"best"` in brackets):

| rows per law | 1% noise | 5% noise | 10% noise |
|---|---|---|---|
| about 100 | 0.20% / 0.18% (0.70%) | 1.3% / 1.2% (3.4%) | 3.4% / 2.6% (6.7%) |
| about 30 | 0.76% / 0.46% (2.0%) | 3.4% / 2.8% (8.4%) | 6.7% / 5.7% (15%) |

- **Test 16's win over PySR came from PySR's default choice.** With `"score"`, PySR finds
  about as many laws as basin at about 100 rows, and more at about 30 rows with
  noise.
- **At about 100 rows the two are level, but on different laws.** Basin keeps every
  product-of-powers law (A) at every noise level, where PySR loses up to 10. PySR finds
  4 to 8 trig and exponential laws (C), which basin cannot write.
- **At about 30 rows PySR is better with noise.** Most of the difference is the 4-input
  laws, which then have 2 levels per input: basin finds 13, 3, 2 and 0 of 27,
  PySR 13, 10, 8 and 7. This is basin's weak spot from Test 16.
- **PySR's answers now predict slightly better than basin's** at every noise level
  and both sizes.
- **`"score"` hardly costs PySR anything.** It loses one law that `"best"` found (I.30.5
  at about 100 rows and 1% noise: it answers `lambd/(d*n)`, the small-angle form of
  `asin(lambd/(n*d))`), and on clean data it finds 26 and 17 more laws.

### Run time against number of rows (Tests 15 and 16)

Every method ran 18 jobs in parallel on the laptop above, so both methods are slower than
on an idle machine. In the Test 16 pilot on an idle machine, PySR took 35–45 s and
basin 3–13 s per run.

| rows per law | basin, median seconds per run | PySR, median seconds per run | PySR slower by | total CPU-hours, basin / PySR |
|---|---|---|---|---|
| about 1000 (Test 15) | 46 | 218 | 4.8× | 34 / 181 |
| about 100 (Test 16) | 17 | 79 | 4.6× | 13 / 55 |
| about 30 (Test 16) | 16 | 61 | 3.8× | 12 / 43 |

- **PySR is slower at every size.** The gap shrinks with fewer rows: about 5× at 1000
  and 100 rows, about 4× at 30.
- **Basin's time barely changes between 100 and 30 rows** (17 against 16 s). Most of
  its time goes on trying combinations of terms and on repeated cross-validation (5 folds
  × 5 repeats), not on the rows themselves.
- **PySR gets faster as rows drop.** It runs a fixed number of search rounds and checks
  every candidate formula on every row, so fewer rows means less work per round.

### Tests 13–14: BACON's own five laws (basin was designed after seeing these)

Exact-law runs. "given" is the authors' own small design and "dense" a larger grid. Same
data for every method, run for run. BACON is the original at its code defaults.

- **PySR 2.7.2 settings:** defaults with `+ - * /` only, maximum size 30, deterministic.
- **Kepler is outside PySR's formulas** with these operators: it cannot write `P^(2/3)`.
- **This is a development check, not a fair contest:** basin was designed knowing
  how BACON and PySINDy failed on these laws.

| law | design | noise | BACON | PySINDy | basin | PySR |
|---|---|---|---|---|---|---|
| ideal gas `V = M(T + 273)/P` | given (27 rows) | none (1 run) | 0 | 0 | **1** | **1** |
|  |  | 1% (20 runs) | 0 | 0 | 0 | **2** |
|  |  | 5% (20 runs) | 0 | 0 | 0 | **3** |
|  |  | 10% (20 runs) | 0 | 0 | 0 | **1** |
|  | dense (216 rows) | none (1 run) | 0 | 0 | **1** | **1** |
|  |  | 1% (20 runs) | 0 | 0 | 0 | **19** |
|  |  | 5% (20 runs) | 0 | 0 | 0 | **5** |
|  |  | 10% (20 runs) | 0 | 0 | 0 | 0 |
| Ohm (wire) `I = T·D²/(2(L + 3))` | given (27 rows) | none (1 run) | **1** | 0 | **1** | 0 |
|  |  | 1% (20 runs) | 0 | 0 | **18** | 0 |
|  |  | 5% (20 runs) | 0 | 0 | **6** | 0 |
|  |  | 10% (20 runs) | 0 | 0 | 0 | 0 |
|  | dense (216 rows) | none (1 run) | **1** | 0 | **1** | 0 |
|  |  | 1% (20 runs) | 0 | 0 | **20** | 6 |
|  |  | 5% (20 runs) | 0 | 0 | **12** | 2 |
|  |  | 10% (20 runs) | 0 | 0 | **3** | 1 |
| Black `T_f = (M₁T₁ + M₂T₂)/(M₁ + M₂)` | given (81 rows) | none (1 run) | **1** | 0 | **1** | **1** |
|  |  | 1% (20 runs) | 0 | 0 | **20** | 8 |
|  |  | 5% (20 runs) | 0 | 0 | 4 | **5** |
|  |  | 10% (20 runs) | 0 | 0 | 0 | 0 |
|  | dense (1296 rows) | none (1 run) | **1** | 0 | **1** | **1** |
|  |  | 1% (20 runs) | 0 | 0 | **20** | 10 |
|  |  | 5% (20 runs) | 0 | 0 | 0 | 0 |
|  |  | 10% (20 runs) | 0 | 0 | 0 | 0 |
| Kepler `D = P^(2/3)` | given (6 rows) | none (1 run) | **1** | 0 | **1** | 0 |
|  |  | 1% (20 runs) | **18** | 0 | 13 | 0 |
|  |  | 5% (20 runs) | **13** | 0 | 6 | 0 |
|  |  | 10% (20 runs) | 5 | 0 | **6** | 0 |
|  | dense (12 rows) | none (1 run) | **1** | 0 | **1** | 0 |
|  |  | 1% (20 runs) | **20** | 0 | 18 | 0 |
|  |  | 5% (20 runs) | 15 | 0 | **18** | 0 |
|  |  | 10% (20 runs) | 4 | 0 | **15** | 0 |
| synthetic `Y = 3/(X + 2)` | given (5 rows) | none (1 run) | **1** | 0 | **1** | **1** |
|  |  | 1% (20 runs) | **8** | 0 | 3 | 0 |
|  |  | 5% (20 runs) | 0 | 0 | 0 | 0 |
|  |  | 10% (20 runs) | 0 | 0 | 0 | 0 |
|  | dense (12 rows) | none (1 run) | **1** | 0 | **1** | **1** |
|  |  | 1% (20 runs) | 2 | 0 | **15** | 0 |
|  |  | 5% (20 runs) | 0 | 0 | **1** | 0 |
|  |  | 10% (20 runs) | 0 | 0 | 0 | 0 |

### Real neuroscience data (exploratory; full write-ups in `findings/`)

- **Three measured datasets,** run on 2026-11-16 with the same method settings as above
  (PySR: the Feynman operators, both `"best"` and `"score"` from one search).
- **Not a pre-registered test.** The rules were set during the session, not written down
  before running.
- **Scoring:** every number is a held-out error. Each method's formula is evaluated on
  cells or neurons it did not fit on, and compared with the measured value.
- **Fits:**
  - Allen: 5 splits of 80/20 over 2,336 cells.
  - NeuroMorpho: 5 fits of 2,000 neurons each, scored on the same 50,000 held-out neurons.
- **basin-free** is basin with `POWER_FREE = True` (see the settings table). It was
  added *after* basin at its defaults failed the wiring test, so
  its wiring result is a development check, like Tests 13–14, not an independent test.

**Held-out median error / RMS log error** (lower is better):

| dataset | basin | basin-free | PySR "best" | PySR "score" | PySINDy | best simple reference |
|---|---|---|---|---|---|---|
| Allen: rheobase ← R, τ, gap | 29.0% / 0.486 | 28.8% / 0.483 | 26.6% / **0.443** | 31.9% / 0.481 | **24.6%** / 0.530 | 25.6% / 0.405 (`515·gap/R`) |
| NeuroMorpho: volume ← S, L, d | **7.8% / 0.232** | **7.8% / 0.232** | 30.1% / 0.531 | 36.9% / 0.500 | 22.9% / 1.371 | 1.5% / 0.290 (`S·d/4`) |
| NeuroMorpho: length ← W, H, D, n | 30.9% / 0.592 | **22.2% / 0.383** | 38.7% / 0.561 | 59.4% / 1.066 | 30.3% / 1.555 | 19.3% / 0.339 (power law, log fit) |

**Readable law found** (fits giving it):

| dataset | basin | basin-free | PySR | PySINDy |
|---|---|---|---|---|
| Allen | `gap/R` 6 of 6 | `gap/R` 3 of 6, `gap^0.8/R` 3 | `gap/R` 0 of 6 | about 300 terms |
| Volume | cylinder law 5 of 5 (`S²/L`) | the same | `c·S·d` 3 of 10 picks | 195 terms |
| Wiring | none (sums of integer terms) | the power law in 5 of 5, plus small extra terms | none | 617 terms |

- **Where basin is clearly ahead: the volume and wiring tests.** Its held-out error is
  far below PySR's and PySINDy's on volume, and below both on wiring with
  `POWER_FREE`. It also finds a readable law in every fit, where PySINDy keeps hundreds of
  terms and PySR changes its answer from fit to fit.
- **Allen: no clear winner on error.**
  - PySINDy has the lowest median error and PySR "best" the lowest log error.
  - Basin is the only method that finds the textbook form `gap/R`, but its fitting
    rule (relative error) pulls the coefficient to 393, against 515 for a log-scale fit.
- **No method beats a simple reference,** except where that reference is the textbook
  law with the wrong constant.
- **`POWER_FREE` is not the default.** It changes the Allen answer to `gap^0.8/R` on
  half the splits. It has not been run on the Feynman benchmark (Tests 15–17), where
  free powers could add wrong candidates.

## Files

| file | what it is |
|---|---|
| `basin/core.py` | the method |
| `basin/__init__.py` | the library API: `discover`, `settings` |
| `hybrid.py` | alias: `import hybrid` gives `basin.core` (used by the runners) |
| `example.py` | the quick-start example: Black's law from a noisy 81-row grid |
| `tests/` | tests (`make test`) |
| `requirements.txt` | the exact package versions used for every result here (Python 3.10) |
| `benchmarks/test15_feynman.csv` | Test 15, every law × method × noise level: runs, exact (fixed and original check), right terms, median prediction error, median seconds |
| `benchmarks/test15_hybrid_pysr_answers.csv` | every basin and PySR answer in Test 15 (4 significant digits), with its verdict |
| `benchmarks/test16_feynman.csv` | Test 16, every grid × law × method × noise level: levels per input, rows, runs, exact (fixed and original check), right terms, median prediction error, median seconds |
| `benchmarks/test16_hybrid_pysr_answers.csv` | every basin and PySR answer in Test 16 (4 significant digits), with its verdict |
| `benchmarks/test17_feynman.csv` | Test 17, every grid × law × noise level for basin, PySR with "best" (Test 16) and PySR with "score": runs, exact (fixed and original check), right terms, median prediction error |
| `benchmarks/test17_pysr_score_answers.csv` | every PySR answer with "score" in Test 17 (4 significant digits), with its verdict |
| `benchmarks/tests12_14_bacon_laws.csv` | Tests 12–14, BACON's five laws, every method |
| `findings/` | write-ups of the real-data runs (Allen, NeuroMorpho volume and wiring, the power switches) |
| `allen/run_allen.py`, `allen/results/` | Allen Cell Types runner, per-fit answers and `allen_scores.csv` |
| `neuromorpho/run_neuromorpho.py`, `neuromorpho/run_wiring.py`, `neuromorpho/results*/` | NeuroMorpho runners (volume, wiring), per-fit answers and scores |
| `data/allen/`, `data/neuromorpho/` | not in the repo (197 MB): the raw downloads (Allen API JSON, NeuroMorpho morphometry pages and merged CSV) |

## Per-law results, Test 15: basin / PySR

Exact-law runs, basin first, then PySR; the larger is in bold. Kind: A one product of
powers, B other algebraic, C with `sin`, `cos`, `exp` or `arcsin`.

| law | kind | formula | no noise (1 run) | 1% (10 runs) | 5% (10 runs) | 10% (10 runs) |
|---|---|---|---|---|---|---|
| I.6.2 | C | `f = exp(-(theta/sigma)**2/2)/(sqrt(2*pi)*sigma)` | 0 / 0 | 0 / 0 | 0 / **1** | 0 / **1** |
| I.6.2a | C | `f = exp(-theta**2/2)/sqrt(2*pi)` | 0 / **1** | 0 / **3** | 0 / **4** | 0 / **1** |
| I.6.2b | C | `f = exp(-((theta-theta1)/sigma)**2/2)/(sqrt(2*pi)*sigma)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.8.14 | B | `d = sqrt((x2-x1)**2+(y2-y1)**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.10.7 | B | `m = m_0/sqrt(1-v**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.12.1 | A | `F = mu*Nn` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.12.2 | A | `F = q1*q2*r/(4*pi*epsilon*r**3)` | 1 / 1 | 10 / 10 | **10** / 9 | 10 / 10 |
| I.12.4 | A | `Ef = q1*r/(4*pi*epsilon*r**3)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.12.5 | A | `F = q2*Ef` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.13.4 | B | `K = 1/2*m*(v**2+u**2+w**2)` | 0 / 0 | 0 / **1** | 0 / 0 | 0 / 0 |
| I.14.3 | A | `U = m*g*z` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.14.4 | A | `U = 1/2*k_spring*x**2` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.15.1 | B | `p = m_0*v/sqrt(1-v**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.15.3t | B | `t1 = (t-u*x/c**2)/sqrt(1-u**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.15.3x | B | `x1 = (x-u*t)/sqrt(1-u**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.16.6 | B | `v1 = (u+v)/(1+u*v/c**2)` | 0 / **1** | 0 / **4** | 0 / **2** | 0 / 0 |
| I.18.4 | B | `r = (m1*r1+m2*r2)/(m1+m2)` | **1** / 0 | **10** / 3 | **10** / 4 | **10** / 3 |
| I.18.12 | C | `tau = r*F*sin(theta)` | 0 / 0 | 0 / **10** | 0 / **10** | 0 / **10** |
| I.18.14 | C | `L = m*r*v*sin(theta)` | 0 / **1** | 0 / **10** | 0 / **10** | 0 / **10** |
| I.24.6 | B | `E_n = 1/2*m*(omega**2+omega_0**2)*1/2*x**2` | 0 / 0 | 0 / 0 | 0 / **3** | 0 / 0 |
| I.25.13 | A | `Volt = q/C` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.26.2 | C | `theta1 = arcsin(n*sin(theta2))` | 0 / **1** | 0 / **10** | 0 / **10** | 0 / **10** |
| I.27.6 | B | `foc = 1/(1/d1+n/d2)` | 1 / 1 | **10** / 9 | 10 / 10 | 10 / 10 |
| I.29.4 | A | `k = omega/c` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.29.16 | C | `x = sqrt(x1**2+x2**2-2*x1*x2*cos(theta1-theta2))` | 0 / 0 | 0 / **1** | 0 / **1** | 0 / 0 |
| I.30.3 | C | `Int = Int_0*sin(n*theta/2)**2/sin(theta/2)**2` | 0 / 0 | 0 / **2** | 0 / **2** | 0 / 0 |
| I.30.5 | C | `theta = arcsin(lambd/(n*d))` | 0 / **1** | 0 / **10** | 0 / **9** | 0 / 0 |
| I.32.5 | A | `Pwr = q**2*a**2/(6*pi*epsilon*c**3)` | 1 / 1 | **10** / 7 | **10** / 8 | **10** / 6 |
| I.34.1 | B | `omega = omega_0/(1-v/c)` | **1** / 0 | 10 / 10 | **10** / 7 | **10** / 7 |
| I.34.8 | A | `omega = q*v*B/p` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.34.14 | B | `omega = (1+v/c)/sqrt(1-v**2/c**2)*omega_0` | 0 / 0 | 0 / **1** | 0 / 0 | 0 / 0 |
| I.34.27 | A | `E_n = (h/(2*pi))*omega` | **1** / 0 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.37.4 | C | `Int = I1+I2+2*sqrt(I1*I2)*cos(delta)` | 0 / 0 | 0 / **2** | 0 / **2** | 0 / 0 |
| I.38.12 | A | `r = 4*pi*epsilon*(h/(2*pi))**2/(m*q**2)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.39.1 | A | `E_n = 3/2*pr*V` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.39.11 | B | `E_n = 1/(gamma-1)*pr*V` | 1 / 1 | 10 / 10 | 10 / 10 | 4 / **10** |
| I.39.22 | A | `pr = n*kb*T/V` | **1** / 0 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.43.16 | A | `v = mu_drift*q*Volt/d` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.43.31 | A | `D = mob*kb*T` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.43.43 | B | `kappa = 1/(gamma-1)*kb*v/A` | 1 / 1 | 10 / 10 | 10 / 10 | 7 / **10** |
| I.47.23 | A | `c = sqrt(gamma*pr/rho)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| I.48.2 | B | `E_n = m*c**2/sqrt(1-v**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| I.50.26 | C | `x = x1*(cos(omega*t)+alpha*cos(omega*t)**2)` | 0 / 0 | 0 / **7** | 0 / **7** | 0 / **8** |
| II.3.24 | A | `flux = Pwr/(4*pi*r**2)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.4.23 | A | `Volt = q/(4*pi*epsilon*r)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.6.11 | C | `Volt = 1/(4*pi*epsilon)*p_d*cos(theta)/r**2` | 0 / 0 | 0 / **10** | 0 / **10** | 0 / **10** |
| II.6.15b | C | `Ef = p_d/(4*pi*epsilon)*3*cos(theta)*sin(theta)/r**3` | 0 / **1** | 0 / **8** | 0 / **9** | 0 / **8** |
| II.8.7 | A | `E_n = 3/5*q**2/(4*pi*epsilon*d)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.8.31 | A | `E_den = epsilon*Ef**2/2` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.10.9 | B | `Ef = sigma_den/epsilon*1/(1+chi)` | 1 / 1 | 10 / 10 | 10 / 10 | 9 / 9 |
| II.11.27 | B | `Pol = n*alpha/(1-(n*alpha/3))*epsilon*Ef` | 0 / 0 | 0 / **1** | 0 / 0 | 0 / 0 |
| II.11.28 | B | `theta = 1+n*alpha/(1-(n*alpha/3))` | **1** / 0 | **10** / 0 | **2** / 0 | 0 / 0 |
| II.13.17 | A | `B = 1/(4*pi*epsilon*c**2)*2*I/r` | 1 / 1 | 10 / 10 | **10** / 9 | 10 / 10 |
| II.13.23 | B | `rho_c = rho_c_0/sqrt(1-v**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| II.13.34 | B | `j = rho_c_0*v/sqrt(1-v**2/c**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| II.15.4 | C | `E_n = -mom*B*cos(theta)` | 0 / 0 | 0 / **4** | 0 / **4** | 0 / **5** |
| II.15.5 | C | `E_n = -p_d*Ef*cos(theta)` | 0 / 0 | 0 / **4** | 0 / **4** | 0 / **5** |
| II.24.17 | B | `k = sqrt(omega**2/c**2-pi**2/d**2)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| II.27.16 | A | `flux = epsilon*c*Ef**2` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.27.18 | A | `E_den = epsilon*Ef**2` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.34.2 | A | `mom = q*v*r/2` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.34.2a | A | `I = q*v/(2*pi*r)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.34.11 | A | `omega = g_*q*B/(2*m)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.34.29a | A | `mom = q*h/(4*pi*m)` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.37.1 | B | `E_n = mom*(1+chi)*B` | 0 / 0 | 0 / **10** | 0 / **10** | 0 / **10** |
| II.38.3 | A | `F = Y*A*x/d` | **1** / 0 | 10 / 10 | 10 / 10 | 10 / 10 |
| II.38.14 | B | `mu_S = Y/(2*(1+sigma))` | 1 / 1 | 10 / 10 | 10 / 10 | **10** / 8 |
| III.4.32 | C | `n = 1/(exp((h/(2*pi))*omega/(kb*T))-1)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| III.4.33 | C | `E_n = (h/(2*pi))*omega/(exp((h/(2*pi))*omega/(kb*T))-1)` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| III.7.38 | A | `omega = 2*mom*B/(h/(2*pi))` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| III.8.54 | C | `prob = sin(E_n*t/(h/(2*pi)))**2` | 0 / 0 | 0 / **2** | 0 / **1** | 0 / 0 |
| III.10.19 | B | `E_n = mom*sqrt(Bx**2+By**2+Bz**2)` | 0 / 0 | 0 / 0 | 0 / **1** | 0 / 0 |
| III.12.43 | A | `L = n*(h/(2*pi))` | **1** / 0 | 10 / 10 | 10 / 10 | 10 / 10 |
| III.13.18 | A | `v = 2*E_n*d**2*k/(h/(2*pi))` | 1 / 1 | 10 / 10 | 10 / 10 | 10 / 10 |
| III.15.12 | C | `E_n = 2*U*(1-cos(k*d))` | 0 / **1** | 0 / **10** | 0 / **9** | 0 / **9** |
| III.15.14 | A | `m = (h/(2*pi))**2/(2*E_n*d**2)` | 1 / 1 | **10** / 8 | **10** / 9 | **10** / 9 |
| III.15.27 | A | `k = 2*pi*alpha/(n*d)` | **1** / 0 | 10 / 10 | 10 / 10 | 10 / 10 |
| III.17.37 | C | `f = beta*(1+alpha*cos(theta))` | 0 / 0 | 0 / **10** | 0 / **10** | 0 / **10** |
| III.21.20 | A | `j = -rho_c_0*q*A_vec/m` | 1 / 1 | **10** / 5 | **10** / 2 | **10** / 2 |
