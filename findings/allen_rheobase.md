# Allen Cell Types: current needed to fire

## Data

- **Source:** Allen Brain Map API, `model::EphysFeature`, all rows (2,336 neurons: 1,920
  mouse, 413 human, 3 unlabelled). Downloaded 2026-11-16 to `data/allen/`. No missing
  values, and all inputs and outputs are positive.
- **Variables:**
  - **I:** `threshold_i_long_square`, rheobase in pA. It is measured with 1 s steps,
    raised in 10 to 20 pA increments.
  - **R:** input resistance in MΩ. The main run uses `ri`, from subthreshold steps. The
    rerun uses `input_resistance_mohm`.
  - **τ:** `tau`, the membrane time constant in ms.
  - **gap:** `threshold_v_long_square − vrest`, in mV.
- **Known answer:** Ohm's law at steady state, I = gap/R. In these units that is
  I = 1000·gap/R. τ should drop out, because a 1 s step reaches steady state.

## Protocol

- Each method fits all cells once (the reported formula), and also 5 times on 80% of the
  cells (5 splits, seed 0).
- The held-out error is the formula evaluated on each split's other 20%, averaged over
  the 5 splits.
- Reference rows: 1000·gap/R as it stands, and c·gap/R with c fitted on each training
  split (least squares in log).

## Results with R = `ri`

| method | formula on all cells | gap/R found (5 splits) | held-out median error | held-out RMS log error | invalid predictions |
|---|---|---|---|---|---|
| known law, c = 1000 | `1000·gap/R` | — | 96% | 0.78 | 0 |
| known law, c fitted | `515·gap/R` | — | **25.6%** | **0.405** | 0 |
| **Hybrid** | **`393·gap/R`** | **5 of 5** | 29.0% | 0.486 | 0 |
| PySR "best" | `6.66e4/(R·√τ)` | 0 of 5 | 26.6% | 0.443 | 0 |
| PySR "score" | `1.90e4/R` | 0 of 5 | 31.9% | 0.481 | 0 |
| PySINDy | about 300 terms | 0 of 5 | 24.6% | 0.530 | 1.8% |

- **The hybrid found the known form in all 6 fits,** with coefficients from 388.6 to
  396.9 and τ dropped. It is the only method that did.
- **The hybrid's coefficient is the optimum of its own fitting rule.** 393.1 is exactly
  the c that minimises mean squared relative error for c·gap/R. Other rules give larger
  values: least squares in log gives 514.6, plain least squares 553.9, and the median
  ratio 509.8.
  - With about 40% scatter, fitting on relative error pulls the coefficient down a lot.
    This is the same bias the main README describes for Test 15, at a much larger noise
    level.
  - This bias accounts for the hybrid's whole gap to the fitted known law (29.0% against
    25.6%).
- **PySR "score" dropped the gap in every fit** (`≈1.9e4/R`).
- **PySR "best" returned a different formula on almost every split:** `R·√τ`, `R + τ`,
  and once `415.9·gap/R + 591.8/τ`.
- **PySINDy has the lowest median error, but it is not the best answer.** It fits
  absolute error, which favours typical-sized cells. It is worse on the log scale, some
  of its predictions are negative currents, and its formula cannot be read.
- **No method clearly beats the fitted known law** on held-out cells. That is consistent
  with gap/R being the real structure and the remaining 25% being scatter.

## The physical coefficient is out of reach

- Every reasonable fit puts c between 390 and 554, against 1000 for Ohm's law: real cells
  need about half the predicted current.
- The likely causes:
  - **The membrane is not linear up to threshold.** Sodium and other channels start
    opening near threshold and add their own current.
  - **Threshold is measured at spike onset,** above the steady voltage the step settles
    at.
- So the textbook law is right about which variables matter and how they combine, and
  wrong about the size of the effect. The 5% coefficient check from the Feynman tests
  fails for every method.

## Rerun with R = `input_resistance_mohm`

| method | answers over the 6 fits | held-out median error |
|---|---|---|
| known law, c fitted | `597·gap/R` | 36.8% |
| Hybrid | `gap/R·(184 + 4481/τ)` (3 splits and all cells), `3047·gap/(R·τ^⅔)` (2 splits) | 26.6% |
| PySR "best" | `2.6e5/(R·τ)` (3 splits), `gap + 2.2e5/(R·τ)` (all cells and 2 splits) | 28.0% |
| PySR "score" | `2.6e5/(R·τ)` | 28.4% |
| PySINDy | about 300 terms | 26.7% |

- **This column fits the law worse.** Its log correlation with 1000·gap/R is 0.73,
  against 0.85 for `ri`, and its free-fit R exponent is −0.80.
- **Every method brings in τ and beats the known law here.** The likely reason is that
  τ = R_true·C carries information about the true resistance that this column measures
  poorly. That points to measurement, not new physics.

## Run time

The 36 jobs (both columns) ran in parallel in 2 min 21 s:

| method | seconds per fit |
|---|---|
| Hybrid | 21 to 28 |
| PySR | about 138 |
| PySINDy | about 7 |

## Reproduce

```bash
P=~/Desktop/bacon-vs-pysindy/.venv-pysr/bin/python     # PySR 2.7.2, pysindy 2.1.0, numpy 2.2.6
cd allen
$P run_allen.py job --method Hybrid --rcol ri --fold all     # methods: Hybrid PySINDy PySR; folds: all 0-4
$P run_allen.py score                                        # writes results/allen_scores.csv
```
