# NeuroMorpho.org: neuron volume from surface, length and diameter

## Data

- **Source:** NeuroMorpho.org API `/api/morphometry`, all 601 pages of 500 (300,495
  neurons, no duplicates). Downloaded 2026-11-16 in 25 minutes, one page at a time.
- **Files:** pages in `data/neuromorpho/pages/`, merged into
  `data/neuromorpho/morphometry.csv`.
- **Quality:** every neuron has a positive surface, volume, length and diameter, so all
  300,495 are used.
  - Some values are extreme and look like tracing or unit errors (a diameter of 1,980 µm,
    a volume of 2×10¹² µm³). They are kept.
- **Variables:**
  - **S:** `surface`, µm²
  - **L:** `length`, total µm
  - **d:** `diameter`, the average, µm
  - **V:** `volume`, µm³
- **Reference law:** the reconstructions are made of cylinders. For one cylinder,
  S = π·d·L and V = (π/4)·d²·L. So V = S·d/4 = (π/4)·d²·L = S²/(4π·L).

## How well the law holds (all 300,495 neurons)

| check | median | middle 50% |
|---|---|---|
| V / (S·d/4) | 1.013 | 1.000 to 1.199 |
| V / ((π/4)·d²·L) | 1.008 | 1.000 to 1.218 |
| S / (π·d·L) | 1.000 | 1.000 to 1.016 |

- **About 48% of neurons match V = S·d/4 exactly.** 61% are within 10%.
- **No neuron falls below the law;** the 1st percentile is exactly 1.0.
  - Volume depends on d², and the average of d² is always at least the square of the
    average d. Neurons with uniform branches sit at 1, and uneven ones sit above (up to
    about 3× at the 99th percentile).
  - So the deviation is real geometry, not noise.
- **The free log-log fit gives V ∝ S^3.15·L^−2.14·d^−1.13.** That looks wrong, but it is
  the law.
  - It equals S²/L × (S/(d·L))^1.15, and S/(d·L) ≈ π is nearly constant.
  - Because S ≈ π·d·L, any amount of that factor fits equally well. Many different-looking
    formulas are the same law, so methods are scored by held-out prediction error rather
    than by matching one form.

## Protocol

- **Splits:** one random permutation (seed 0) of all neurons. The first 50,000 form the
  held-out test set.
- **Fits:** 5, each on its own 2,000 neurons (disjoint, none from the test set).
- **Scoring:** every formula is evaluated on the 50,000 test neurons and compared with
  their recorded V. Results are averaged over the 5 fits.
- **Reference rows:** S·d/4 as it stands, and c·S·d with c fitted on each training
  sample (least squares in log).

## Results

| method | typical formula | cylinder law found | constants | median error | within 10% | RMS log error | invalid predictions |
|---|---|---|---|---|---|---|---|
| reference S·d/4 | `0.25·S·d` | — | 1 | **1.5%** | 62% | 0.290 | 0% |
| reference, c fitted | `0.287·S·d` | — | 1 | 14.9% | 29% | 0.256 | 0% |
| **Hybrid** | **`0.0861·S²/L`** | **5 of 5** | 1 | 7.8% | **79%** | **0.232** | 0% |
| PySR "best" | e.g. `S·(d/(3.475 − cos(L)) + 105.4·d/L)` | 0 of 5 | 2 to 4 | 30.1% | 27% | 0.531 | 0.03% |
| PySR "score" | `0.26`–`0.34·S·d` | 3 of 5 | 1 to 4 | 36.9% | 28% | 0.500 | 0.03% |
| PySINDy | 195 terms | 0 of 5 | 195 | 22.9% | 29% | 1.371 | 13% |

The hybrid's answers, fit by fit:

| fit | formula |
|---|---|
| 0 | `0.08606·S²/L` |
| 1 | `0.08571·S²/L` |
| 2 | `0.08618·S²/L` |
| 3 | `0.08606·S²/L` |
| 4 | `0.02705·S³/(L²·d)` |

- **Every hybrid fit is the cylinder law.** Fit 4 equals about 0.085·S²/L, because
  S/(L·d) ≈ π.
- **Its coefficient is 8% above the textbook 1/(4π) = 0.0796.** No neuron falls below the
  law, so fitting on relative error is pulled upward by the neurons above it.
- **On median error the textbook law wins** (1.5% against 7.8%), because about half the
  neurons match it exactly. No fitted coefficient can hit that cluster and the tail at
  once. The fitted-coefficient reference shows the same trade-off.
- **On "within 10%" and RMS log error, the hybrid beats the textbook law.** The reason is
  in [neuromorpho_cylinder_forms.md](neuromorpho_cylinder_forms.md).

PySR's answers, fit by fit:

| fit | "score" | "best" |
|---|---|---|
| 0 | `0.318·S·d` | `0.08741·S²/(L − √S)` |
| 1 | `0.261·S·d` | `107.7·S·(0.0007538·S + d)/L` |
| 2 | `0.007053·S·(cos(2016.7·d) + 4.672)/(L/S − 0.1806/d)` | same |
| 3 | `0.503·S` (diameter dropped) | `d²·(√L·d·(−d·sin(L·d − d) + d + 1.626) + 0.7714·L)` |
| 4 | `0.340·S·d` | `S·(d/(3.475 − cos(L)) + 105.4·d/L)` |

- **"score" found c·S·d in 3 fits,** with coefficients up to 36% off 0.25.
- **"best" never found the law,** and three of its answers contain `sin` or `cos` of a
  length or diameter.
- **PySINDy returned 195 terms in every fit,** and 13% of its held-out predictions are
  negative or infinite.

## Compared with Allen

This data is much cleaner: about 8% scatter, against 25% on Allen. Here the hybrid is
clearly ahead of PySR and PySINDy on held-out error as well as on finding the form. On
Allen it was ahead only on the form.

## Run time

The 15 jobs ran in parallel in 84 s:

| method | seconds per fit |
|---|---|
| Hybrid | about 19 |
| PySR | about 73 |
| PySINDy | about 3 |

## Reproduce

```bash
P=~/Desktop/bacon-vs-pysindy/.venv-pysr/bin/python
cd neuromorpho
$P run_neuromorpho.py job --method Hybrid --fit 0     # methods: Hybrid PySINDy PySR; fits 0-4
$P run_neuromorpho.py score                           # writes results/neuromorpho_scores.csv
```
