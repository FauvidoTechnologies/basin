# NeuroMorpho: S²/L beats S·d, and the hybrid chose it

## The question

- **Three forms are equal for a cylinder:** V = S·d/4 = S²/(4π·L) = (π/4)·d²·L.
- **They differ on real neurons,** because a neuron is many cylinders of different
  diameters.
- **Which form predicts best?** Each form is scored with the same two coefficients, so
  that only the form differs:
  - the textbook coefficient;
  - the coefficient the hybrid's own rule (least squared relative error) gives, averaged
    over the 5 training samples.
- **Test set:** the same 50,000 held-out neurons as in
  [neuromorpho_volume.md](neuromorpho_volume.md). 47.8% of them match S·d/4 exactly and
  47.9% match S²/(4πL) exactly.

## Results

On the 26,088 held-out neurons that are not exact cylinders:

| form | coefficient | median error | within 10% | RMS log error |
|---|---|---|---|---|
| S·d | textbook 0.25 | 16.1% | 27% | 0.401 |
| **S²/L** | **textbook 1/(4π)** | **14.5%** | **32%** | **0.360** |
| d²·L | textbook π/4 | 17.4% | 25% | 0.455 |
| S·d | fitted 0.2716 | 9.1% | 53% | 0.353 |
| **S²/L** | **fitted 0.08599** | **7.8%** | **58%** | **0.315** |
| d²·L | fitted 0.8552 | 10.4% | 49% | 0.407 |

On all 50,000 held-out neurons:

| form | coefficient | median error | within 10% | RMS log error |
|---|---|---|---|---|
| S·d | textbook | 1.5% | 62% | 0.290 |
| S²/L | textbook | 1.3% | 64% | 0.260 |
| d²·L | textbook | 1.4% | 61% | 0.329 |
| S·d | fitted | 8.6% | 76% | 0.261 |
| S²/L | fitted | 8.1% | 78% | 0.234 |
| d²·L | fitted | 8.9% | 73% | 0.300 |

## What this shows

- **S²/L wins on every measure,** with both coefficients and on both sets of neurons. The
  margin is moderate, not dramatic.
- **The ranking follows how much each form uses the diameter.** S²/L does not use it at
  all, S·d uses it once, and d²·L uses it twice (as d²).
  - That suggests NeuroMorpho's "average diameter" is a worse summary of a neuron's
    thickness than its surface and length are.
  - A possible reason: d is averaged over tracing points or compartments, while volume
    sums d² × length branch by branch.
  - This is not verified against NeuroMorpho's documentation.
- **The hybrid wrote the best form.** It returned `0.0861·S²/L` in 4 of 5 fits, matching
  the best c·S²/L above (0.08599). Its fifth answer, `0.02705·S³/(L²·d)`, is the same
  form given S ≈ π·d·L.
  - Its step 4 compares candidate laws by cross-validated relative error, and S²/L
    scored best.
  - So it never wrote S·d itself. That is the form a person would write, but not the one
    that predicts best here.
- **PySR "score" was the method that returned S·d** (3 of 5 fits, coefficients 0.26 to
  0.34). PySINDy's term list contains S·d, S²/L and d²·L, but it never kept one of them
  alone.

## Open questions

- **How NeuroMorpho computes "diameter".** Is it averaged per point, per compartment, or
  length-weighted? This would confirm or rule out the explanation above.
- **Whether S²/L stays better within each archive or tracing system.** NeuroMorpho
  combines many labs, which may compute diameters differently.
