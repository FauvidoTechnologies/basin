# Hybrid: BACON proposes terms, regression decides

A law finder for tables of measurements: given inputs `x₁ … xₙ` and one measured output
`y`, it returns a formula `y = N(x)/D(x)`. It combines ideas from BACON (Langley's
discovery system, as rebuilt by Miller & Banerjee 2025) with SINDy-style sparse
regression (PySINDy), and decides between candidate laws by cross-validation.

It was built and benchmarked in
[bacon-vs-pysindy](https://github.com/pUrGe12/bacon-vs-pysindy) (Tests 13–15). This
directory holds a standalone copy of it, so it can be developed and tested against PySR
on its own.

## Where this copy came from

- `hybrid.py` is copied from `bacon-vs-pysindy/hybrid.py` at commit `a53d1ce`. That file
  (SHA-256 `b418a97a…`) is the version recorded in all 79 hybrid run records of Test 15,
  unchanged since it was frozen for Test 13.
- **One change:** the original imported a list of 25 STLSQ thresholds
  (`np.logspace(-3, 3, 25)`) from `sindy_runner.py`. The list is now written into
  `hybrid.py`, so the file stands alone. Nothing else differs (`diff` shows only those
  lines).
- **Checked:** the copy and the original were run on the same data for four Test 15
  setups, and gave identical formulas, equal to the saved Test 15 answers:
  - I.12.1 at 10% noise: `0.9751*Nn*mu`
  - I.18.4 at 5% noise: `(0.9975*m1*r1 + 0.9982*m2*r2)/(1.0*m1 + 1.009*m2)`
  - I.47.23 with no noise: `1.0*gamma**0.5*pr**0.5/rho**0.5`
  - II.11.28 at 1% noise: `(0.6733*alpha*n + 1.0)/(-0.3292*alpha*n + 1.0)`
- Any change from here on should be a named switch, with the benchmark version as the
  default, so the results below stay reproducible.

## Use

```bash
pip install -r requirements.txt      # Python 3.10; the versions used for every result below
python example.py                    # Black's law from 81 noisy rows, about 10 s
```

```python
import numpy as np, sympy as sp
from hybrid import discover

# X: one row per measurement, one column per input; y: the measured output
law = discover(X, y, [sp.Symbol(v, positive=True) for v in ["M1", "T1", "M2", "T2"]])
law, info = discover(X, y, symbols, details=True)   # also: proposed powers, CV error
```

## How it works

Given the table (every step uses all rows):

1. **Propose powers, BACON's way.** For each input `x_j`, group the rows where every
   other input has the same value (BACON's "slices" of a factorial experiment) and fit
   `log|y| = (one intercept per slice) + p·log x_j` on all slices at once. Every simple
   fraction (denominator 1, 2 or 3, |p| ≤ 3, not 0) within 3 standard errors of `p`
   (at least 0.03; at most the 2 nearest) becomes a candidate term `x_j^p`. So does the
   product of each input's nearest fraction (BACON's invariant, e.g. `T·D²`). Nothing is
   accepted or rejected here. If the rows have no slices, one fit of `log|y|` on all
   `log x_j` is used instead.
2. **Term list:** `1`, every input, every product of two inputs (squares included), plus
   the proposed terms.
3. **Fractions as a linear problem (as in SINDy-PI):** `y = N/D` with `N` and `D` sparse
   sums of terms, fitted as `y·D = N` with one `D` coefficient fixed to 1. Candidate laws:
   every law with up to 4 terms in total (screened by a quick fit, best 30 per size),
   plus larger ones (up to 5 free coefficients) from PySINDy's STLSQ and from forward
   selection.
4. **Decide by held-out error.** Each candidate is scored by 5-fold cross-validation
   repeated over 5 random splits, on mean squared relative error, with coefficients
   refined by least squares on the relative error. Pick the fewest terms whose error is
   within one standard error of the best; ties go to lower total power.

| setting | value | meaning |
|---|---|---|
| `MAX_FREE` | 5 | most free coefficients in a law (N and D together) |
| `EXHAUSTIVE_TERMS` | 4 | every law with up to this many terms is tried |
| `SCREEN_KEEP` | 30 | per law size, candidates kept after the quick screen |
| `POWER_SE`, `POWER_FLOOR` | 3, 0.03 | step 1 tolerance: 3 standard errors, at least 0.03 |
| `POWER_KEEP` | 2 | at most 2 proposed fractions per input |
| `FOLDS`, `REPEATS`, `CV_SEED` | 5, 5, 0 | cross-validation (leave-one-out below 10 rows) |
| `CV_FLOOR` | 1e-10 | errors closer than this count as equal |

**What it cannot write:** `sin`, `cos`, `exp`, `log`; powers with a denominator above 3 or
above 3 in size; laws with more than 5 free coefficients; products of more than two
inputs unless step 1 proposes them; fractional powers of inputs that are not all
positive. Step 1 is skipped if any `y` is exactly 0.

## Benchmarks

All results are from bacon-vs-pysindy. The rules for each test were written down before
it ran, and the hybrid was not changed after seeing any result. "Exact law" means the
right terms with every number within 5% of its true value (equivalent rearrangements
count). Noise is Gaussian on the measured output with sd = level × |y|.

### Test 15: the AI Feynman benchmark (79 laws, the hybrid was not designed on these)

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
| **Hybrid** | 43 of 79 | 43 of 79 | 42 of 79 | 41 of 79 |
| **PySR** | 42 of 79 | 50 of 79 | 50 of 79 | 49 of 79 |

**The same, split by kind of law.** The split was made after the run, from the true
formulas:
- **A:** one product of powers, e.g. `q1*q2/(4*pi*epsilon*r**2)`.
- **B:** other algebraic laws: sums, fractions with sums, roots of sums.
- **C:** laws with `sin`, `cos`, `exp` or `arcsin`.

| kind of law | noise | BACON (original) | BACON-robust | PySINDy | Hybrid | PySR |
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

**Hybrid against PySR, law by law.** Each count is the number of laws where that method
finds the exact law in more runs than the other. The sign test is two-sided. The "all 79"
rows were the test's decision rule. The splits by kind were made after the run, so they
describe the results rather than test them.

| laws | noise | Hybrid finds the law in more runs | PySR finds it in more runs | same | sign test |
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

| laws | noise | runs | within 0.1%: Hybrid | PySR | within 1%: Hybrid | PySR | within 5%: Hybrid | PySR |
|---|---|---|---|---|---|---|---|---|
| all 79 | no noise | 79 each | 43 | 66 | 53 | 73 | 61 | 75 |
|  | 1% noise | 790 each | 423 | 483 | 528 | 718 | 618 | 766 |
|  | 5% noise | 790 each | 0 | 261 | 483 | 554 | 580 | 749 |
|  | 10% noise | 790 each | 0 | 214 | 0 | 470 | 592 | 675 |
| 59 algebraic | no noise | 59 each | 43 | 53 | 50 | 57 | 55 | 57 |
|  | 1% noise | 590 each | 423 | 377 | 500 | 557 | 558 | 580 |
|  | 5% noise | 590 each | 0 | 206 | 472 | 429 | 540 | 570 |
|  | 10% noise | 590 each | 0 | 176 | 0 | 380 | 543 | 523 |

**The hybrid's numbers shrink with noise.** On runs that find the exact law, the median
prediction error is:

| noise | Hybrid: exact runs | median error | PySR: exact runs | median error |
|---|---|---|---|---|
| no noise | 43 | 0.00% | 42 | 0.00% |
| 1% noise | 430 | 0.04% | 522 | 0.02% |
| 5% noise | 422 | 0.57% | 517 | 0.09% |
| 10% noise | 410 | 2.27% | 491 | 0.21% |

On simple products the hybrid answers `0.9927*Nn*mu` at 5% noise and `0.9751*Nn*mu` at
10% (`F = mu*Nn`; the same on `m*g*z`, `q/C`, `Ef*q2`), in every run. PySR answers `Nn*mu`.

- **The cause is the hybrid's fitting rule.** Step 4 fits on error relative to the noisy
  `y`, `(prediction − y)/y`. Dividing by a noisy `y` pulls every coefficient down by about
  2σ²: 0.5% at 5% noise, 2% at 10%.
- **The verdicts don't change.** These answers are still within the 5% check, so this
  does not change any count above. But at around 15–20% noise the bias alone would push
  coefficients past 5%.
- **Part of PySR's edge here is the benchmark.** Many Feynman laws have the coefficient
  1, which PySR can write with no fitted number at all.

**Run time** (all methods running 18 jobs in parallel on a 22-thread laptop at about
1.2 GHz. Rerun alone on this laptop, four hybrid runs took 3 to 9 times less time):

| inputs | laws | Hybrid, median seconds per run | PySR, median seconds per run |
|---|---|---|---|
| 1 | 1 | 6 | 131 |
| 2 | 15 | 17 | 102 |
| 3 | 36 | 43 | 362 |
| 4 | 27 | 69 | 149 |
| all | 79 | 46 (total 34 CPU-hours) | 218 (total 181 CPU-hours) |

### Run time against number of rows (Tests 15 and 16)

Test 16 reruns Test 15 on smaller grids: the same 79 laws, settings, noise levels and
seeds, with about 100 or about 30 rows per law (fewer levels per input). Its law counts
will be added here when it finishes. Times are from finished runs, with every method
running 18 jobs in parallel on the laptop above, so both methods are slower than on an
idle machine. In the Test 16 pilot on an idle machine, PySR took 35–45 s and the hybrid
3–13 s per run.

| rows per law | laws timed | Hybrid, median seconds per run | PySR, median seconds per run | PySR slower by | total CPU-hours, Hybrid / PySR |
|---|---|---|---|---|---|
| about 1000 (Test 15) | 79 of 79 | 46 | 218 | 4.8× | 34 / 181 |
| about 100 (Test 16) | 79 of 79 | 17 | 79 | 4.6× | 13 / 55 |
| about 30 (Test 16, partial) | Hybrid 63, PySR 53 of 79 | 18 | 59 | 3.3× | not final |

*The about-30 row is partial, as of 10 November 2026, 23:09.*

- **PySR is slower at every size.** The gap shrinks with fewer rows: about 5× at 1000
  and 100 rows, about 3× at 30.
- **The hybrid's time barely changes between 100 and 30 rows** (17 against 18 s). Most of
  its time goes on trying combinations of terms and on repeated cross-validation (5 folds
  × 5 repeats), not on the rows themselves.
- **PySR gets faster as rows drop.** It runs a fixed number of search rounds and checks
  every candidate formula on every row, so fewer rows means less work per round.

### Tests 13–14: BACON's own five laws (the hybrid was designed after seeing these)

Exact-law runs. "given" is the authors' own small design and "dense" a larger grid. Same
data for every method, run for run. BACON is the original at its code defaults.

- **PySR 2.7.2 settings:** defaults with `+ - * /` only, maximum size 30, deterministic.
- **Kepler is outside PySR's formulas** with these operators: it cannot write `P^(2/3)`.
- **This is a development check, not a fair contest:** the hybrid was designed knowing
  how BACON and PySINDy failed on these laws.

| law | design | noise | BACON | PySINDy | Hybrid | PySR |
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

### What the benchmarks say about the hybrid against PySR

- **Algebraic laws: level.** On the 59 laws without `sin`, `cos`, `exp` or `arcsin`, the
  two find about the same laws at every noise level (after-the-run split, not
  significant).
  - The hybrid is more reliable on some laws: 10 of 10 runs at every noise level, where
    PySR gets 3–4 of 10 on I.18.4 `(m1*r1 + m2*r2)/(m1 + m2)` and 2–5 on III.21.20
    `-rho_c_0*q*A_vec/m`.
  - PySR finds II.37.1 `mom*(1 + chi)*B` in 10 of 10 runs at every noise level, and
    I.16.6 on clean data. The hybrid finds neither in any run. `chi` gets no proposed
    power, so the term `mom*chi*B` is not in its list.
- **Trig and exponential laws: PySR only.** The hybrid cannot write them. This is the
  whole of PySR's lead in the overall counts.
- **Noise: the hybrid's law count barely moves** (43, 43, 42, 41 laws). But its
  coefficients shrink with noise (above), so PySR's answers are closer to the true law
  at 5–10% noise.
- **Clean data: PySR sometimes adds a tiny extra term**, e.g.
  `omega_0*v/(c - v + 7.18e-8) + omega_0`. That costs it 11 laws that it finds at 1%
  noise.
- **Small designs (5–27 rows, Tests 13–14): the hybrid holds up better.** PySR builds
  long formulas that fit the few points. Neither finds the ideal gas law, because
  `T + 273` and its first-order twin fit equally well over the narrow range of `T`.
- **Cost: the hybrid needs about a fifth of PySR's CPU time** (34 against 181 CPU-hours
  for Test 15, 13 against 55 at about 100 rows). The gap narrows to about a third at 30
  rows (partial Test 16 timings, above).

## Ideas for testing the hybrid against PySR

As in bacon-vs-pysindy: one test at a time, with rules written down before running and
both methods frozen. Rough cost: about 4 hours per 1,000 PySR runs on this laptop with 18
parallel jobs, about 1 hour per 1,000 hybrid runs.

**No new runs: re-score the saved Test 15 answers**

1. **Extrapolation.** Score every saved answer at points outside the training ranges
   (e.g. each input up to twice its range). An exact law keeps working outside the range.
   A formula that only fits inside it stops working, e.g. PySR's
   `1.012*m_0*(1 + exp(-c/v))` for I.10.7. This shows how much each method's "close but
   not exact" answers are worth. It costs minutes.

**New runs, both methods unchanged**

2. **Random samples instead of a grid.** Test 15 gave every method a factorial grid,
   which is BACON's kind of data and the hybrid's home turf. Standard symbolic-regression
   benchmarks (SRBench, SRSD) sample inputs at random, so there are no slices, and step 1
   falls back to one global fit. This asks whether the hybrid's strength depends on
   designed experiments. Same 79 laws and row counts.
3. **Noise of constant size** (sd = level × RMS of y, the SRBench convention) instead of
   noise proportional to |y|. The current noise suits the hybrid's log fits and relative
   errors. With constant noise, small values of `y` get large relative noise and can
   change sign. Pairs naturally with idea 2 as one test with four conditions: grid or
   random, proportional or constant noise.
4. **Higher noise: 15%, 20%, 30%.** This tests the bias above. The prediction is that
   the hybrid starts missing laws it has the right terms for, at around 15–20%.
5. **Fewer rows: 20, 50, 100 instead of about 1000.** Real experiments are often this
   small. Tests 13–14 suggest PySR overfits small tables and the hybrid does not; this
   would test that on laws neither was designed on.
6. **Constants that are not 1.** Many Feynman laws can be written with no fitted number,
   which PySR does exactly. Two ways to remove that edge:
   - Multiply each law by a random constant.
   - Use SRSD-Feynman (Matsubara et al.), which keeps the physical constants and more
     realistic ranges.

   This also tests the hybrid with very small or large values, e.g. ε₀ ≈ 8.85e-12.
7. **Irrelevant inputs.** Add 1–3 inputs that do not affect `y`; SRSD has "dummy
   variable" versions for this. Both methods must leave them out. The hybrid's term list
   grows with the square of the number of inputs.
8. **More inputs: the 21 Feynman laws with 5–9 inputs**, held out so far.
   - 8 of them use trig or exp: I.12.11, I.40.1, I.41.16, II.11.17, II.35.18, II.35.21,
     III.9.52, III.14.14.
   - This tests how the hybrid's search scales, and whether its 5-coefficient limit
     blocks laws.
9. **PySR given more:**
   - 10 times the iterations: does more search find more laws, or add more extra
     terms?
   - Its "accuracy" model selection instead of "best", which would show whether the
     clean-data extra terms come from model selection.
   - Only the hybrid's operators on the algebraic laws.

   This checks that PySR is not handicapped by its default settings.

**Changes to the hybrid (each a switch, frozen before the test)**

10. **Remove the coefficient bias:** refine on error relative to the fitted value instead
    of the noisy `y`, repeating the fit. Test together with idea 4.
11. **`sin`, `cos`, `exp` shapes in step 1.** Propose `g(x_j)` for a small set of shapes,
    the same way as powers. This could reach 9 of the 20 trig laws: I.18.12, I.18.14,
    II.6.11, II.6.15b, II.15.4, II.15.5, III.17.37, I.37.4, I.6.2a. The other 11 need the
    inside of the function found too (e.g. `ω·t` or `θ/σ`).
12. **Ablations: switch off one part at a time** — the BACON step, the fraction form, or
    cross-validation (in-sample error instead). This shows which part makes the
    difference against PySR before building more on it (e.g. ABACUS ideas).
13. **The BACON step inside PySR:** give PySR the proposed terms `x_j^p` as extra inputs.
    This shows whether BACON's proposals help a genetic search too.

**Suggested order:**
1. Idea 1: it is free.
2. Ideas 2 + 3 as the next test, run on fresh noise seeds. It would also repeat the
   after-the-run "level on algebraic laws" result as a proper test.
3. Ideas 4 + 10, then idea 11.

## Files

| file | what it is |
|---|---|
| `hybrid.py` | the method (copy of bacon-vs-pysindy's, see above) |
| `example.py` | small example: Black's law from a noisy 81-row grid |
| `requirements.txt` | package versions used for every result above |
| `benchmarks/test15_feynman.csv` | Test 15, every law × method × noise level: runs, exact (fixed and original check), right terms, median prediction error, median seconds |
| `benchmarks/test15_hybrid_pysr_answers.csv` | every hybrid and PySR answer in Test 15 (4 significant digits), with its verdict |
| `benchmarks/tests12_14_bacon_laws.csv` | Tests 12–14, BACON's five laws, every method |

## Per-law results, Test 15: hybrid / PySR

Exact-law runs, hybrid first, then PySR; the larger is in bold. Kind: A one product of
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
