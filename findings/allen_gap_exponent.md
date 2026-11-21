# Allen: the gap exponent is 0.79 to 0.86, not 1

## What was fitted

- **The free log-log fit** assumes I = k·gap^a·R^b·τ^c, so
  log I = log k + a·log gap + b·log R + c·log τ, solved by least squares on all 2,336
  cells (R = `ri`).
- **With τ:** a = 0.79, b = −0.99, c = −0.27.
- **Without τ:** a = 0.86, b = −1.12.
- **What a means:** at equal R and τ, a cell with twice the gap needs 2^0.79 ≈ 1.7× the
  current, not 2×.

## Why it should be 1

At steady state a passive membrane is a resistor, so I = ΔV/R. Current is proportional
to voltage, which means exponent 1 on the gap and −1 on R.

## Two explanations that fit equally well

1. **Real physics.** Near threshold, sodium and other channels open and help the
   depolarisation along. A cell with a larger gap crosses more of that region, so it
   needs proportionally less extra current. That would give an exponent below 1.
2. **Measurement error in the gap.** Noise in an input pulls its fitted exponent toward 0,
   by about var(true)/(var(true) + var(noise)).
   - The gap has a narrow spread: most cells are between 21.4 and 44.4 mV, and the sd of
     log gap is 0.236.
   - So a few mV of error is enough:

| measurement noise on the gap | a true exponent of 1 would come out as |
|---|---|
| 1 mV | about 0.98 |
| 2 mV | about 0.92 |
| 3 mV | about 0.82 |
| 4 mV | about 0.68 |

   - R spans 65 to 340 MΩ (sd of log R 0.502), so similar noise barely moves its
     exponent. That is consistent with b landing at −0.99 to −1.12.

## The data cannot separate them

- **The reverse fit gives an upper bound of 4.6.** Regressing log gap on log I and log R
  and inverting gives 4.6 as the bound on the gap exponent. With errors in the
  variables, the true value lies between the forward fit (0.86) and that bound.
- **The range includes 1 by a wide margin.**
- **What would settle it:** repeated measurements on the same cells, to estimate the
  gap's measurement noise. The Allen feature table does not have them.

**Conclusion:** the gap exponent is somewhere around 0.8 to 1. That is consistent both
with Ohm's law plus measurement noise and with a real nonlinearity near threshold.
