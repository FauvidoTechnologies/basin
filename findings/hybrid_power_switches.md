# Hybrid side runs: a larger fraction set, and free powers

## Why

- **On the wiring test,** the hybrid proposed no powers. Its estimates (W^0.263,
  H^0.573, D^0.168, n^0.416, each ±0.01) are not within max(3 SE, 0.03) of any k/d with
  d ≤ 3.
- **So it lost to a plain power law fitted in log space** (see
  [neuromorpho_wiring.md](neuromorpho_wiring.md)).
- **Two named switches were added to `hybrid.py`.** The defaults are the benchmark
  version: `example.py` gives byte-identical output, and a rerun of wiring fit 0
  reproduces its saved answer exactly.

| switch | default | side run | effect |
|---|---|---|---|
| `POWER_DENOMINATORS` | `(1, 2, 3)` | `(1, 2, 3, 4, 5, 6)`, called **Hybrid-den6** | step 1 also tries quarters, fifths and sixths |
| `POWER_FREE` | `False` | `True`, called **Hybrid-free** | when no fraction is close enough, step 1 proposes the estimate itself, rounded to `POWER_FREE_DIGITS = 2` decimals |

- **Nothing else changed.** Splits, held-out sets and scoring are the same as in the
  original tests.
- **Rerun with the switches:** wiring, NeuroMorpho volume, and Allen with R = `ri`.

## Wiring (the target): both switches help, and free powers come close to the best power law

| method | median error | RMS log error | within 25% |
|---|---|---|---|
| free power law, log fit (reference) | **19.3%** | **0.339** | **61%** |
| **Hybrid-free** | **22.2%** | **0.383** | **56%** |
| Hybrid-den6 | 25.7% | 0.432 | 49% |
| Cuntz, c fitted (reference) | 29.4% | 0.508 | 43% |
| Hybrid (benchmark) | 30.9% | 0.592 | 40% |
| PySINDy | 30.3% | 1.555 | 44% |
| PySR "best" | 38.7% | 0.561 | 34% |

The answers:

| fit | Hybrid-den6 | Hybrid-free |
|---|---|---|
| 0 | `2.42·D^(1/6)·H^(3/5)·W^(1/4)·n^(2/5) / (1 + 1.3e-5·same)` | `2.37·D^0.17·H^0.57·W^0.26·n^0.42 + 15.9` |
| 1 | `2.49·D^(1/6)·H^(3/5)·W^(1/4)·n^(2/5) / (1 + 3.0e-5·same)` | `2.65·D^0.16·H^0.57·W^0.26·n^0.42 / (1 + 2.4e-5·same)` |
| 2 | `(3.12·D^(1/6)·H^(3/5)·n^(2/5) + 2.03·W)/(1 − 0.17·D^(1/6))` | `2.61·D^0.14·H^0.56·W^0.29·n^0.42` |
| 3 | `(4.80·D + 0.66·H·n)/n^(2/5)` | `2.67·D^0.16·H^0.56·W^0.29·n^0.40 / (1 + 0.0029·H^0.56)` |
| 4 | `2.58·D^(1/6)·H^(3/5)·W^(1/4)·n^(2/5) / (1 + 0.0030·H^(3/5))` | `2.23·D^0.16·H^0.56·W^0.27·n^0.43 + 1.19·W` |

- **Both switches let step 1 propose every input's power.** BACON's invariant (the
  product of all four) then becomes a candidate term, and it appears in 9 of the 10
  answers.
- **Most answers are that power law with a small correction:** a constant, an extra W,
  or a mild saturation 1/(1 + ε·…).
- **Hybrid-free fit 2 is exactly the plain power law,** with exponents within 0.02 of the
  log fit's.
- **The remaining gap to the log-fit reference (22.2% against 19.3%) is probably the
  fitting rule.** The hybrid fits relative error, which with this much scatter pulls the
  coefficients away from the log-scale optimum, as on Allen. The small extra terms may
  also overfit. Neither cause is verified.
- **Cost:** about 90 s per fit, against 68 s for the benchmark hybrid.

## NeuroMorpho volume: no harm

| method | median error | RMS log error | within 10% |
|---|---|---|---|
| Hybrid (benchmark) | 7.81% | 0.2316 | 78.7% |
| Hybrid-den6 | 8.05% | 0.2340 | 78.2% |
| Hybrid-free | 7.81% | 0.2316 | 78.7% |

- **Hybrid-free gave the same 5 answers as the benchmark.**
- **Hybrid-den6 gave the same answers on fits 0 to 3.** On fit 4 it wrote `0.0859·S²/L`
  in place of the equivalent `0.0271·S³/(L²·d)`.
- **Both are 2 to 3 times slower** (40 to 70 s per fit, against 13 to 30 s), because the
  term list is longer.

## Allen (R = `ri`): same accuracy, but the textbook form is lost on most splits

| split | Hybrid (benchmark) | Hybrid-den6 | Hybrid-free |
|---|---|---|---|
| all cells | `393.1·gap/R` | `393.1·gap/R` | `393.1·gap/R` |
| 0 | `396.9·gap/R` | `940.3·gap^(3/4)/R` | `396.9·gap/R` |
| 1 | `396.5·gap/R` | `704.2·gap^(5/6)/R` | `712.3·gap^0.83/R` |
| 2 | `390.8·gap/R` | `924.1·gap^(3/4)/R` | `390.8·gap/R` |
| 3 | `392.7·gap/R` | `931.3·gap^(3/4)/R` | `2758·gap^(2/3)/(R·τ^0.28)` |
| 4 | `388.6·gap/R` | `774.3·gap^(4/5)/R` | `774.3·gap^0.8/R` |

| method | held-out median error | RMS log error |
|---|---|---|
| Hybrid (benchmark) | 29.0% | 0.486 |
| Hybrid-den6 | 29.7% | 0.492 |
| Hybrid-free | 28.8% | 0.483 |

- **The finer vocabulary lets the hybrid write the gap exponent the data shows** (0.75
  to 0.83, close to the free log fit's 0.79 to 0.86), instead of Ohm's 1.
  - Hybrid-den6 does so on 4 of 5 splits, Hybrid-free on 3 of 5.
  - On all cells, both still choose `gap/R`.
- **Predictions barely change** (within 0.7 points of median error). Cross-validation
  cannot tell gap^0.8 from gap^1 on this data, so which one is chosen depends on the
  split.
- **Whether that is a step back depends on the goal.** See
  [allen_gap_exponent.md](allen_gap_exponent.md): the 0.8 may be measurement error
  pulling the exponent down rather than physics. The benchmark hybrid's coarse fractions
  act as a prior toward simple laws, which here gives back the textbook law. A finer set
  instead reports the exponent the noisy data shows.

## Summary

| test | Hybrid-den6 | Hybrid-free |
|---|---|---|
| Wiring | better: 30.9% to 25.7% | best: 30.9% to 22.2%, close to the log-fit power law (19.3%) |
| Volume | the same (one equivalent rewrite) | identical |
| Allen | gap/R kept on 1 of 5 splits | gap/R kept on 2 of 5 splits |
| Cost | 1.3 to 3 times slower | 1.3 to 3 times slower |

- **`POWER_FREE` is the better switch:** more accurate on wiring, fewer changes on Allen,
  and no change on volume.
- **Neither is the default.** All Feynman, BACON and Allen results above it were produced
  with the defaults. Adopting a switch as the default would need a rerun of the Feynman
  benchmarks (Tests 15 to 17), where a denser vocabulary could add wrong proposals.

## Reproduce

```bash
P=~/Desktop/bacon-vs-pysindy/.venv-pysr/bin/python
cd neuromorpho && $P run_wiring.py job --method Hybrid-free --fit 0      # or Hybrid-den6
$P run_neuromorpho.py job --method Hybrid-free --fit 0
cd ../allen && $P run_allen.py job --method Hybrid-free --rcol ri --fold 0
```
