# Busy beaver growth census

An exploratory study, run on 2026-11-23: what basin's fitted laws reveal about small
Turing machines. Each machine is simulated, its tape width is recorded over time, and
basin is asked for a law `span(t)`. The laws sort the machines into a few growth
classes. Their constants turn out to identify the mechanism each machine uses.

**Not a pre-registered test.** The classification rules (step 5 of the pipeline) were
fixed while looking at partial results, after the first version misclassified some
machines. Every number below is recomputed from the committed files by
`python busy_beaver/analyze.py`.

**Not a halting test either.** A fitted growth law says how a machine behaves in the
window it was fitted on, not what happens after. The BB(5) champion fits `1.92·√t`
cleanly and halts at step 47,176,870 (see [the champion](#the-bb5-champion-c-derived-from-its-rule)).

## Findings in short

1. **A few sharp growth classes, in almost the same shares for 3, 4 and 5 states.**
   About 83% linear, 13% bounded, 2–4% √t, 0.4% logarithmic, plus a rare t^(1/3) class
   from 4 states up. Nothing in between in the samples, but see the
   [limits](#limits-and-open-questions).
2. **The constants are exact numbers that identify the mechanism:**
   - **Linear speeds** are small fractions (1/3, 1, 1/5, 1/7, …). A parity rule explains
     why odd denominators dominate.
   - **Bouncers:** c² in `span ≈ c·√t` is a simple fraction for 60% of machines, against
     12% for random values over the same range.
   - **Counters** count in exact algebraic bases. Of the 297,225 counters in the BB(5)
     database, 15,054 count in **base φ (Fibonacci / Zeckendorf counters)**. Others count
     in the plastic number, its square, the tribonacci constant, Narayana's cows constant
     and the silver ratio. All of these are **Pisot numbers**.
3. **The BB(5) champion's constant follows from its Collatz-like rule.** span/√t swings
   between exactly 4/√5 and 2 once per factor of 25/9 in time. The average, 1.9236,
   matches basin's 1.920.
4. **How badly basin fits detects hard machines, but not how they end.** Every BB(6)
   holdout fits more roughly than 99% of ordinary machines (AUC 1.000). Within the
   holdouts, the cryptids are not special, and probviously-halting cryptids look like
   probviously non-halting ones (p = 0.49).

---

## Data

### Machines

| set | where from | size | used |
|---|---|---|---|
| BB(5) | bbchallenge seed database, `all_5_states_undecided_machines_with_global_header.zip` from `https://docs.bbchallenge.org/` (254 MB zipped, 2.66 GB unzipped) | 88,664,064 records: 14,322,029 time-limit, then 74,342,035 space-limit | uniform sample of 20,000 (seed 0), 0.023%; the counter census uses all 14,322,029 time-limit records |
| BB(4) | own tree enumeration, `small/enum.c` | 436,278 non-halting machines | uniform sample of 20,000 (seed 0), 4.6% |
| BB(3) | own tree enumeration | 2,678 non-halting machines | all, 100% |
| BB(6) holdouts | `https://wiki.bbchallenge.org/w/images/5/5e/BB6_holdouts_608.txt` (shared by mxdys, 30 Oct 2026), copied to `bb6/BB6_holdouts_608.txt` | 608 (up to equivalence) | all |
| BB(6) cryptids | bbchallenge wiki, Cryptid and BB(6) pages, read 2026-11-23 | 20 cryptids and 6 potential cryptids | all; listed in `results/named_machines.csv` |
| BB(6) champion | bbchallenge wiki | `1RB1RA_1RC1RZ_1LD0RF_1RA0LE_0LD1RC_1RA0RE` | runs about 10↑↑10↑↑10↑↑8 steps |
| BB(5) top halters | bbchallenge wiki, BB(5) page, Top Halters table | 20 machines, 2,358,064 to 47,176,870 steps | all; each halts at exactly its published step count in our simulator |
| random 6-state | `bb6/rand6.c`: A0 = 1RB, one halting slot, other transitions uniform | 20,000 distinct machines that visit all 6 states and survive 2^16 steps (2.6% of draws) | 5,000 fitted |
| BB(5) halters | `halt/enum_halt.c`: the same tree, cut off at 2^16 steps | 118 halters that halt between 2^12 and 2^16 steps | all |

- **Record format (BB(5) database):**
  - Header: 30 bytes, starting with three big-endian uint32 counts.
  - Records: 30 bytes each, 5 states × 2 read symbols × (write, move 0=R 1=L, next
    state 1–5, 0 = undefined/halt).
  - All 88,664,064 machines are proven never to halt (BB(5), 2024).
- **Tree enumeration (3 and 4 states):**
  - Same tree normal form as bbchallenge: first transition `1RB`, new states
    introduced in order.
  - A machine still running after S(n)+1 steps (22 and 108) can never reach an undefined
    transition, so it never halts.
  - Check: the deepest halting leaves are exactly 6, 21 and 107 steps for n = 2, 3, 4.
  - Like bbchallenge's tree, it never produces a machine with every transition defined.
    A 3-state machine with all 6 transitions set appears only in the 4-state tree,
    with state D unused.
- **Counts are per record.** Records that differ only in transitions the machine never
  uses are counted separately, as in the bbchallenge database.

### Measurement

- **Dependent variable:** `span`, the cells from the leftmost to the rightmost cell
  visited so far, inclusive.
- **Independent variable:** `t`, the step count. 49 checkpoints at
  `t_k = round(2^(10 + k/4))`, k = 0…48 (2^10 to 2^22 = 4,194,304 steps). The fits use
  k = 8…48 (41 points, 2^12 to 2^22).
- **Simulator:** `census/sim.c` reads database records; `bb6/simtxt.c` reads machine
  strings and handles any number of states. Both reproduce the BB(5) champion's
  47,176,870 steps exactly. None of the non-halting machines halts within the run.

## Pipeline

1. **Simulate** each machine from a blank tape, state A, for 2^22 steps. Record the span
   at the 49 checkpoints.
2. **Bounded:** if the span is the same at all 41 fit checkpoints, the class is
   bounded and basin is not called.
3. **Fit twice** with `basin.discover(X, y, [name], details=True)` at default settings
   (`MAX_FREE=5`, `EXHAUSTIVE_TERMS=4`, 5×5 cross-validation, power denominators 1, 2, 3):
   - span against `t`;
   - span against `u = log₂ t`. basin cannot write `log`, so `a·log t + b` becomes
     `a·u + b` here.
4. **Roughness** is basin's `cv_rel_mse`, the cross-validated relative error of the
   chosen law.
5. **Choose a law** (rules set during the session):
   - A `u` law is only allowed if the span grew slower than t^0.4 between 2^20 and 2^22.
     Without this rule, 3 clean √t machines went to log, because a 5-coefficient
     rational in `u` beat `c·√t` on error.
   - If the two errors are within a factor of 2, keep the law with fewer coefficients.
     Without this rule, 116 √t machines went to log.
   - Otherwise keep the lower error.
6. **Class:** a `u` law gives `log`. Otherwise use d ln span / d ln t of the law at
   t = 2^22, rounded to 0, 1/3, 1/2, 2/3 or 1 if it is within 0.04, and `other` if not.
7. **Rate:**
   - speed `f(t)/t` for linear laws;
   - `f(t)/t^p` for t^p laws;
   - `f(22.5) − f(21.5)`, cells per doubling of t, for log laws;
   - the span itself for bounded machines.

## Results

### Growth classes (window 2^12 to 2^22 steps)

| growth class | 3 states (2,678) | 4 states (20,000) | 5 states (20,000) | 5 states, time limit (3,259) | 5 states, space limit (16,741) |
|---|---|---|---|---|---|
| linear (translated cyclers) | 84.2% | 83.2% | 83.1% | – | 99.3% |
| bounded (cyclers) | 13.5% | 13.4% | 12.8% | 78.4% | – |
| √t (bouncers) | 2.0% | 3.0% | 3.7% | 19.6% | 0.7% |
| log (counters) | 0.4% | 0.4% | 0.3% | 1.9% | – |
| t^(1/3) | – | 0.02% (4) | 0.02% (4) | 0.1% | – |
| other | – | – | 1 machine | – | 1 machine |

- **95% margins** with 20,000 machines: about ±0.5 points on linear and bounded, ±0.26
  on √t, ±0.08 on log.
- **The rare classes are poorly pinned down.** 4 machines in 20,000 puts t^(1/3)
  anywhere from 0.005% to 0.05%.
- **The time-limit / space-limit split is a split by class.** That follows from how the
  database was built: a linear machine reaches 12,289 cells long before 47M steps.
- **The one "other" machine** is `1RB---_1LC1RA_0RB0RD_1RC1LE_0LD1RC`. Its span
  multiplies by a repeating cycle of factors (×1.97, 1.46, 1.42, 1.42) every 32× in
  time, which averages t^0.51. It is √t with a log-periodic wobble that basin has no
  term for.
- Results per machine: `results/bb{3,4,5}_census.csv`.

### Constants

**Linear speeds.**
- A translated cycler repeats a p-step pattern that shifts d cells, so its speed is
  d/p.
- The head moves every step, so **p ≡ d (mod 2)**. Speed 1/3 can be 1 cell per 3
  steps, but 1/2 needs 2 cells per 4 steps, a longer and rarer pattern.
- BB(5) sample counts: 1/3: 2,758, 1: 2,665, 1/5: 2,064, 1/7: 1,266, 1/4: 965,
  1/2: 856, 1/9: 847, 1/11: 583, 1/6: 564, 3/11: 527.
- Speeds a/b with a+b even (the period can be b): 11,465 machines; with a+b odd (the
  period must be at least 2b): 3,362.

**Bouncers, `span ≈ c·√t`.**
- If the region grows by a cells per round trip and the head sweeps at h cells per
  step, then t ≈ a·n²/h after n round trips and span ≈ a·n. So **c² = a·h**.
- In the BB(5) sample, c² falls within 0.2% of a fraction with denominator ≤ 6 for
  60.5% of the 749 √t machines. Random values spread over the same range do that
  11.5% of the time. At denominator ≤ 4 and 0.1% tolerance it is 28.3% against 2.9%.
- The commonest value is c² = 2 (c = √2): one new cell at each end per round trip, at
  full speed. Next come 1, 4/3, 4, 3/2, 8/3, 2/3, 3.

**t^(1/3) machines.**
- Three 4-state machines have law coefficient c = 1.4424 ≈ 3^(1/3).
- Two 5-state machines share `2.2895·t^(1/3) − 0.979`, and 2.2895³ = 12.00.
- A third 5-state machine has `2.892·t^(1/3) − 2.36`; 2.892³ = 24.2, which is not clean.
- Hypothesis: each round of growth costs time ∝ span², so c³ = 3 × (a simple ratio).

**Counters: exact bases.** These come from a separate census, not the 20,000 sample:

1. Scan all 14,322,029 time-limit records for 2^16 steps (`base/scan.c`).
2. Treat a slope of span between 2^14 and 2^16 in (0, 0.25] as a counter: 297,225
   records. (11,229,321 records have constant span; 2.77M have slope 0.4–0.6.)
3. Run each counter to 2^24 steps (`base/incrbin.c`) and record every step at which
   the span grows.
4. Measure the time factor per new cell over the last 6k increments after 2^14 steps.
   Whole multiples of 6 cancel repeating patterns of length 1, 2 and 3.

A counter in base β with k cells per digit gains a cell each time t grows by
β^(1/k).

| time factor per new cell | what it is | records | distinct growth sequences |
|---|---|---|---|
| √2 | binary, 2 cells per digit | 140,307 | 31,636 |
| 2 | binary, 1 cell per digit | 81,492 | 11,011 |
| 2^(1/3) | binary, 3 cells per digit | 31,439 | 12,770 |
| **φ = 1.618034** | **Fibonacci (Zeckendorf), x² = x + 1** | **15,054** | **4,235** |
| 2^(1/4) | binary, 4 cells per digit | 6,880 | 3,883 |
| √3 | ternary, 2 cells per digit | 5,631 | 2,024 |
| 2^(2/3) | binary, 3 cells per 2 digits | 3,102 | 1,332 |
| 3^(1/3) | ternary, 3 cells per digit | 1,623 | 851 |
| √φ | Fibonacci, 2 cells per digit | 614 | 388 |
| ψ = 1.465571 | Narayana's cows constant, x³ = x² + 1 | 476 | 309 |
| ρ = 1.324718 | plastic number (Padovan), x³ = x + 1 | 338 | 209 |
| √(1+√2) | silver ratio (Pell), 2 cells per digit | 302 | 192 |
| ρ² = 1.754878 | plastic number squared | 184 | 134 |
| 2^(1/6), 3^(1/4), 2^(1/5) | more cells per digit | 138, 124, 12 | |
| tribonacci = 1.839287 | x³ = x² + x + 1 | 101 | 53 |
| φ^(1/3) | Fibonacci, 3 cells per digit | 27 | 23 |
| unidentified | | 9,381 (3.2%) | 6,333 |

- **How tight the matches are:** median relative error 7×10⁻⁵ for φ (15,195 records
  within 10⁻³), 2×10⁻⁵ for tribonacci, 1.2×10⁻⁵ for ρ², 2.4×10⁻⁵ for Narayana,
  3.2×10⁻⁵ for silver, and 2.1×10⁻⁴ for ρ, which converges slowly.
- **In the 20,000 sample**, run to 2^32 steps, the four φ machines match to 10⁻⁶
  (`results/bb5_sample_counters_increments_2e32.txt`).
- **"Distinct growth sequences"** means distinct exact lists of increment times. Equal
  lists are almost certainly the same behaviour, so this is a ceiling on the number of
  distinct behaviours.
- **What base φ means.**
  - Positions are worth Fibonacci numbers 1, 2, 3, 5, 8, 13, …, each about φ times the
    previous.
  - No two adjacent digits are 1, because "011" equals "100". Counting 1 to 8:
    1, 10, 100, 101, 1000, 1001, 1010, 10000.
  - Counting to N takes about log_φ N digits. Time grows with the count, so a cell is
    added each time t grows by φ.
- **Pisot hypothesis.** Every identified base is a Pisot number: an algebraic integer
  above 1 whose other roots are all smaller than 1 in absolute value. The plastic
  number is the smallest one. Frougny showed that carries in base β can be done by a
  finite automaton when β is a Pisot number. A small Turing machine is a finite
  automaton sweeping a tape. So small machines may only be able to count in Pisot
  bases. This is **untested**: a counter in a non-Pisot base would refute it.
- **Was this known?** The bbchallenge wiki's Counter page mentions only "binary
  counters, ternary counters, and so on". We found no page on Fibonacci, Zeckendorf or
  non-integer-base counters (checked 2026-11-23: Counter, Collatz-like, Zoology,
  5-state winner; `/wiki/Fibonacci` is a 404). Discord was not checked.
- **Unidentified:** the largest clusters are at 1.361 (4,646 records) and 1.848 (189),
  which match nothing in the candidate list. Some slightly-off values near 1.414
  (about 1,200 records) look like √2 counters that have not settled by 2^24.
- Results: `results/bb5_counter_bases.csv.gz` (`db_index` is the record's index in the
  database).

### The BB(5) champion: c derived from its rule

The rule, from the bbchallenge wiki and Pascal Michel:

- C(n) = `0^∞ <A 1^n 0^∞`, and the blank tape is C(0).
- C(3k) → C(5k+6) in 5k² + 19k + 15 steps.
- C(3k+1) → C(5k+9) in 5k² + 25k + 27 steps.
- C(3k+2) → halt.
- Trajectory: 0, 6, 16, 34, …, 7366, 12284, then halt.

`champ/trace.c` confirms that at every predicted time the tape has exactly the form
C(n) with the published n, that the span is n+1, and that the machine halts at
47,176,870.

The derivation:

- Each step multiplies n by about 5/3 and costs about (5/9)·n² steps. Summing the
  geometric series gives t ≈ (5/16)·n² at each configuration.
- **So at configurations span/√t → 4/√5 = 1.7889.** Measured: 1.7891.
- **Within a step**, the wiki's sub-rule C(m,n) → C(m+5, n−1) in 2m+8 steps gives
  m ≈ √(5τ) and span ≈ N + (2/5)·m, after time τ into a step that starts at C(N). So
  span/√t = (1 + (2/√5)·s) / √(5/16 + s²), with s = √τ / N.
- **Its maximum is exactly 2**, at s = √5/8, which is 22% of the way through each step
  in log t. Measured: peaks of 2.0003 to 2.003, at 21–25%.
- **One cycle is a factor of 25/9 = 2.78 in time.** Measured: ×2.78 to ×2.80.
- **The average over log t is 1.9236.** basin fitted **1.920**.
- **This wobble is all of the champion's "roughness"** (cv 1.5×10⁻³, rougher than 99%
  of non-halting √t machines). It comes from macro steps that multiply the size by a
  fixed factor, which makes growth repeat on a log scale. Bouncers add a fixed amount
  per sweep and grow smoothly.

Data: `results/champion_trace.txt`.

### Roughness: hard machines, BB(6), halters

Roughness is compared within each growth class. Fits: `results/bb6_fits.csv`.

- **Holdouts vs ordinary machines.** 100% of the 608 BB(6) holdouts are rougher than
  the 99th percentile of 5,000 random 6-state machines (AUC 1.000) and of the BB(5)
  census (AUC 0.999).
- **The holdouts' class mix:** no simple law 37.2%, √t 32.4%, log 16.9%, t^(1/3) 8.6%,
  nearly flat 3.1%, t^(2/3) 1.5%, linear 0.3%. Random 6-state machines are 86.9%
  linear and 11.9% bounded.
- **Holdout class labels are less reliable.** The share with a clean law is 75% for
  √t, 44% for t^(1/3) and 16% for log. The 9 apparent t^(2/3) machines are curving
  rationals, not clean powers.
- **Cryptids:** every cryptid in the √t class (including Antihydra, `1.38·√t + 5.7`)
  is rougher than all random √t machines and more than 98.5% of BB(5) non-halting √t
  machines.
- **Within the holdouts, cryptids are spread from the 20th to the 98th percentile.**
  Antihydra is at 36%, the BB(6) champion at 53%.
- **Halting vs non-halting cryptids:** holdout percentile medians 71 (n=10, probviously
  halting) against 55 (n=9, probviously non-halting). Mann-Whitney p = 0.49: no
  difference.
- **BB(5) top halters:** the 19 in the √t class have c = 1.37–1.92 (mostly one family
  near 1.9, halting around 11.8M steps). All are rougher than 98.7% of BB(5) non-halting
  √t machines.
- **Typical BB(5) halters.** The 118 halters that halt between 2^12 and 2^16 steps
  were fitted on 2^10–2^12, alongside 3,000 non-halters on the same window.
  - Halters: √t 52.5%, other 28%, t^(1/3) 8.5%, log 8.5%, t^(2/3) 2.5%, and no linear
    or bounded ones (those loop forever).
  - Within √t they are about 8× rougher than non-halters (median cv 7.7×10⁻⁴ against
    9.1×10⁻⁵, p = 3×10⁻¹²).
  - Data: `results/bb5_halters_2e12_to_2e16.txt`,
    `results/halters_vs_nonhalters_window_2e10_2e12.jsonl`.

The conclusion: roughness separates hard from easy machines. It does not separate
cryptids from other hard machines, or halting from non-halting.

## Limits and open questions

- **Gaps between classes.** No BB(5) machine in the sample has a slope between 0.15 and
  0.25 or between 0.55 and 0.95. With 20,000 machines that only bounds those gaps below
  about 0.015% (rule of three), which could still be about 13,000 records. The rare
  t^(1/3) class shows such things exist.
  - Planned test: raw slopes (C simulator only, no basin) for all 436,278 BB(4)
    machines and 1,000,000 BB(5) machines. Estimated about 20 minutes on 20 cores.
- **Window.** Every class describes 2^12 to 2^22 steps only.
  - Planned check: rerun the same machines for longer (e.g. to 2^26) and see whether
    the classes and constants hold. This is most important for the holdouts' t^(1/3),
    t^(2/3) and "other" labels.
- **Pisot hypothesis:** look for a counter in a non-Pisot base, and identify the 1.361
  and 1.848 clusters.
- **Bouncers:** derive c² = a·h from the transition table for a few machines, and
  check it against the fitted c.
- **Cryptids:** derive Antihydra's c from its rule, as was done for the champion.
- **Baseline:** the random 6-state baseline is uniform over transition tables, not
  over the TNF tree, so its class mix is not comparable to the BB(3)–BB(5) shares.
- **Share stability across n** (83/13/3%) may partly be because bigger trees contain
  smaller machines with unused states.

## Reproducing

Run from inside `busy_beaver/`, with the repository root on `PYTHONPATH` so that
`import basin` works. Compile each `.c` file with `gcc -O2 -o <name> <name>.c -lm`.
Run times are on 20 cores. Large intermediate files are git-ignored (`*.out`, `*.bin`,
`bbdata/`).

| step | command (from the folder named) | time |
|---|---|---|
| download | `curl -o bbdata/db.zip https://docs.bbchallenge.org/all_5_states_undecided_machines_with_global_header.zip` | – |
| BB(5) sample | `census/`: `python sample.py ../bbdata/db.zip 20000 sample` | 5 s |
| simulate | `census/`: `./sim sample.bin <start> <end> 22 part.out` per slice, concatenated to `sim22.out` | 30 s |
| fit | `census/`: `python fit_all.py` (writes `fits.jsonl`) | 75 min |
| classify | `census/`: `python classify.py` (writes `classified.json`) | 1 min |
| BB(3), BB(4) | `small/`: `./enum 3 21 bb3.bin`, `./enum 4 107 bb4.bin`, sample 20,000 of BB(4) (seed 0), simulate with `../census/sim`, `python fit_small.py <sim> <out> <workers>`, then `python ../census/classify.py <sim> <fits> <bin> <idx.npy> <out> <n>` | 90 min |
| BB(6) | `bb6/`: list machine strings (`special.csv`, `BB6_holdouts_608.txt` last field, `./rand6 6 20000 1 65536 rand6.txt`), `./simtxt <list> 22 <out>`, `python fit_any.py 20 special:special22.out:all holdout:holdouts22.out:all random6:rand6_22.out:5000`, `python classify6.py` | 45 min |
| halters | `halt/`: `./enum_halt 5 65536 4096 20 <p> h5_<p>.txt` for p = 0…19, simulate with `../bb6/simtxt <list> 16`, `python fit_window.py 18 0 8 halt_fits.jsonl halter:<sim>:25:all nonhalt:../census/sim22.out:49:3000` | 60 min |
| counters | `base/`: extract the first 14,322,029 records to `bbdata/timelimit.bin`, `./scan` per slice, select slope (0, 0.25], `./incrbin ../bbdata/timelimit.bin cand.int64 <a> <b> 24 <out>` per slice | 35 min |
| champion | `champ/`: `./trace > trace.txt` | 1 s |
| page | `census/`: `python build_page.py` (writes `atlas.html`) | – |

`analyze.py` reads only `results/` and prints every number quoted above.

## Files

| path | what it is |
|---|---|
| `analyze.py` | recomputes every number in this file from `results/` |
| `census/sim.c` | simulator for database records: span and ones at log-spaced checkpoints |
| `census/sample.py` | uniform sample of records, streamed from the zip |
| `census/fit_all.py` | basin fits on `t` and `log₂ t` (BB(5) sample) |
| `census/classify.py` | law choice, growth class and rate; takes paths and n for BB(3) and BB(4) |
| `census/build_page.py`, `template.html`, `class_text.json`, `page_meta.json` | the atlas page |
| `small/enum.c` | tree enumeration of non-halting n-state machines (checked: S(2), S(3), S(4)) |
| `small/fit_small.py` | the fits, skipping the `log₂ t` fit when the slope is above 0.4 |
| `bb6/simtxt.c` | simulator for machine strings, any number of states |
| `bb6/rand6.c` | random machines with one halting slot that visit every state |
| `bb6/fit_any.py`, `bb6/classify6.py` | fits and classes for the named machines, holdouts and random machines |
| `bb6/special.csv`, `bb6/BB6_holdouts_608.txt` | named machines and the holdout list |
| `halt/enum_halt.c` | halting leaves of the tree up to a step cutoff (checked against `enum.c`, and split across processes) |
| `halt/fit_window.py` | fits on a fixed checkpoint window |
| `base/scan.c`, `base/incrbin.c`, `base/incr.c` | counter census: span scan, and span-increment times |
| `champ/trace.c` | the champion's configurations and span/√t |
| `results/bb{3,4,5}_census.csv` | per machine: class, law, cv, rate, span |
| `results/bb6_fits.csv` | BB(6) named machines, holdouts and random machines; also the BB(5) top halters |
| `results/bb5_counter_bases.csv.gz` | all 297,225 BB(5) counters with machine string, time factor per new cell and base |
| `results/champion_trace.txt` | configurations and span at 2,000 log-spaced times |
| `results/named_machines.csv` | cryptids, potential cryptids, BB(6) champion, BB(5) top halters, with sources' statuses |
| `results/bb5_halters_2e12_to_2e16.txt`, `results/halters_vs_nonhalters_window_2e10_2e12.jsonl` | the halter comparison |
| `results/bb5_sample_counters_increments_2e32.txt` | span-increment times to 2^32 steps for the sample's 63 log-class machines |

The interactive atlas built from the first census (classes, histograms and space-time
diagrams) was published as a private claude.ai artifact; `census/build_page.py`
regenerates it.

## Sources

- bbchallenge seed database: https://docs.bbchallenge.org/
- Cryptids: https://wiki.bbchallenge.org/wiki/Cryptid
- BB(6), holdouts and champion: https://wiki.bbchallenge.org/wiki/BB(6),
  https://wiki.bbchallenge.org/wiki/Holdouts_lists
- BB(5) top halters: https://wiki.bbchallenge.org/wiki/BB(5)
- Champion's rule: https://wiki.bbchallenge.org/wiki/5-state_busy_beaver_winner,
  https://bbchallenge.org/~pascal.michel/beh#tm52a
- Counters: https://wiki.bbchallenge.org/wiki/Counter
- Frougny, C. (1992), "Representations of numbers and finite automata", *Mathematical
  Systems Theory* 25, 37–60.
