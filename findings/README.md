# Findings on real data

Hybrid vs PySINDy vs PySR on measured datasets rather than synthetic laws. The method
settings are the benchmark ones from bacon-vs-pysindy, unchanged, with no tuning for
these datasets:

- **Hybrid:** `hybrid.py` with its benchmark defaults.
- **PySINDy:** STLSQ with the Feynman term list (integer powers up to ±3, total degree ≤ 5).
- **PySR:** 2.7.2 with the Feynman operators, maxsize 40, deterministic. One search per
  fit, with both picks recorded: `"best"` (PySR's default) and `"score"`.

Every error is held-out. Each method's formula is evaluated on cells or neurons it did
not fit on, and compared with their measured values.

| file | dataset | question |
|---|---|---|
| [allen_rheobase.md](allen_rheobase.md) | Allen Cell Types, 2,336 neurons | current needed to fire ← input resistance, time constant, voltage gap |
| [allen_gap_exponent.md](allen_gap_exponent.md) | same | why the free fit's gap exponent is 0.79 to 0.86, not 1 |
| [neuromorpho_volume.md](neuromorpho_volume.md) | NeuroMorpho.org, 300,495 neurons | volume ← surface, length, diameter |
| [neuromorpho_cylinder_forms.md](neuromorpho_cylinder_forms.md) | same | S²/L beats S·d, and why the hybrid wrote S²/L |
| [neuromorpho_wiring.md](neuromorpho_wiring.md) | same, 221,899 3D neurons | dendrite length ← width, height, depth, branch points (Cuntz 2012) |
| [hybrid_power_switches.md](hybrid_power_switches.md) | all three | side runs: fractions up to sixths (`POWER_DENOMINATORS`), and free powers (`POWER_FREE`) |

## Short version

| dataset | Hybrid | PySR | PySINDy |
|---|---|---|---|
| Allen, R = `ri` | `393·gap/R` in 6 of 6 fits, the known form | never gap/R (`1.9e4/R`, `R·√τ`, …) | about 300 terms |
| NeuroMorpho volume | `0.086·S²/L` (cylinder law) in 5 of 5 fits, 7.8% held-out median error | c·S·d in 3 of 10 picks, 30–37% median error | 195 terms, 13% invalid predictions |
| NeuroMorpho wiring | no power law: exponents not simple fractions, so none proposed. Sums like `(0.19·H·n + 3.9·W)/(1 + 0.0001·n²)`, 30.9% median error | `√`-products or a single input, 39–59% | 617 terms, 30.3%, 12% invalid |

On the wiring test **all three lose to a plain power law fitted in log space** (19.3%
median error).

Code: also `neuromorpho/run_wiring.py`, results in `neuromorpho/results_wiring/`.

Code: `allen/run_allen.py` and `neuromorpho/run_neuromorpho.py`. Results: `allen/results/`
and `neuromorpho/results/`. Raw data: `data/allen/` and `data/neuromorpho/`, downloaded
2026-11-16.
