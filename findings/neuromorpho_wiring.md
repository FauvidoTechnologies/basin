# NeuroMorpho: dendrite length from bounding box and branch points

## Question and data

- **Question:** total dendrite length L ← width W, height H, depth D (µm) and number of
  branch points n.
- **Literature:** Cuntz, Mathy & Häusser (PNAS 2012) derive L = (3/(4π))^(1/3)·V^(1/3)·n^(2/3)
  from optimal wiring, where V is the spanning volume.
- **Here, V is the bounding box W·H·D,** a crude stand-in. W, H and D are given as separate
  inputs, so the methods have to find any combination themselves.
- **Neurons:** the NeuroMorpho morphometry table (see
  [neuromorpho_volume.md](neuromorpho_volume.md)), filtered to depth ≥ 5 µm (3D
  tracings) and n ≥ 1. That leaves 221,899 neurons.
- **Protocol:** as in the volume test. 50,000 neurons are held out; there are 5 fits on
  2,000 neurons each (seed 0); every formula is evaluated on the held-out neurons.
- **Reference rows:**
  - Cuntz's formula with its textbook constant (0.6204)
  - Cuntz's formula with its constant fitted
  - a free power law c·W^a·H^b·D^c·n^e fitted by least squares in log, on each training
    sample

## Results (50,000 held-out neurons, mean of 5 fits)

| method | typical formula | constants | median error | within 10% | within 25% | RMS log error | invalid predictions |
|---|---|---|---|---|---|---|---|
| **free power law (log fit)** | **`3.0·W^0.27·H^0.57·D^0.16·n^0.42`** | 5 | **19.3%** | **28%** | **61%** | **0.339** | 0% |
| Cuntz, c fitted | `1.73·(W·H·D)^(1/3)·n^(2/3)` | 1 | 29.4% | 17% | 43% | 0.508 | 0% |
| Cuntz, textbook c | `0.620·(W·H·D)^(1/3)·n^(2/3)` | 1 | 63.2% | 2% | 4% | 1.144 | 0% |
| Hybrid | `(0.19·H·n + 3.9·W)/(1 + 0.0001·n²)` | 3 to 4 | 30.9% | 15% | 40% | 0.592 | 0% |
| PySR "best" | e.g. `6.23·H^0.25·n^0.5·(D + 2W)^0.5` | 1 to 2 | 38.7% | 13% | 34% | 0.561 | 0% |
| PySR "score" | `24.6·W`, `17.7·D`, `W·√n` | 0 to 1 | 59.4% | 6% | 17% | 1.066 | 0% |
| PySINDy | 617 terms (the whole term list) | 617 | 30.3% | 20% | 44% | 1.555 | 11.6% |

## Findings

- **No method found a power law. All three lose clearly to a 5-constant power law fitted
  in log space,** which has about 2/3 of their median error and the lowest log error by
  far. This is the first real-data test in this folder where the hybrid does not find the
  best simple law.
- **Why the hybrid missed it:** its step 1 proposes powers only near simple fractions.
  - On fit 0 its own log fit gave W^0.263, H^0.573, D^0.168, n^0.416, each with a
    standard error of about 0.01.
  - None is within max(3 SE, 0.03) of a fraction with denominator 1, 2 or 3. For example,
    0.416 is 0.08 from 1/3 and from 1/2, and 0.168 is 0.17 from 1/3.
  - So no powers were proposed (`proposed_powers = {}`), and the term list held only 15
    integer terms.
  - The hybrid then did the best it could with sums and ratios of those:
    `3.49·D + 0.048·H·n + 2.16·W` on fit 0, and `(a·H·n + b·W)/(1 + c·n²)` on fits 1 to 4.
  - This is the "What it cannot write" limit from the main README, met on real data.
    Here the true exponents are not simple fractions at all.
- **PySR "best"** found square-root products, such as `6.23·H^0.25·n^0.5·(D + 2W)^0.5`,
  which are rough approximations of the power law.
- **PySR "score"** reduced the answer to one or two inputs (`24.6·W`, `17.7·D`). That is
  the worst of all methods apart from the textbook constant.
- **PySINDy** kept its whole 617-term list. It has the second-best median error but the
  worst log error, and 11.6% of its predictions are invalid.

## What the data says about the Cuntz law

- **The size exponents sum to 1, but that is a units check, not evidence for Cuntz.** The
  free fit's W, H and D exponents sum to 0.26 + 0.57 + 0.16 ≈ 0.99, the same as the
  combined fit's V^0.34. L is a length and n has no units, so any law with consistent
  units must have the length exponents summing to 1. Cuntz's V^(1/3) satisfies this, but
  so would any other consistent law.
- **The specific predictions do not hold on bounding boxes:**
  - **The axes are not equal** (0.27, 0.57, 0.16, against 1/3 each). Depth matters least,
    probably because the z axis of a tissue slice is thin, compressed and cut. Height
    matters most, probably because many neurons are traced with their long axis along y.
    Neither explanation is verified.
  - **The branch exponent is about 0.42, not 2/3.** The paper's n is the number of
    branch points, used here as `n_bifs`. A crude box in place of spanning volume, and a
    mix of cell types, may account for some of this.
- **This is not a test of Cuntz's theory as published,** which uses spanning volume and
  was checked within cell types. It shows only that the law does not carry over to
  NeuroMorpho's bounding boxes across all cells.

## Run time

The 15 jobs ran in parallel in 2 min:

| method | seconds per fit |
|---|---|
| Hybrid | about 68 |
| PySR | about 113 |
| PySINDy | about 45 |

## Reproduce

```bash
P=~/Desktop/bacon-vs-pysindy/.venv-pysr/bin/python
cd neuromorpho
$P run_wiring.py job --method Hybrid --fit 0     # methods: Hybrid PySINDy PySR; fits 0-4
$P run_wiring.py score                           # writes results_wiring/wiring_scores.csv
```
