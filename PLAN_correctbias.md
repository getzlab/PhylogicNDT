# Plan: Make CorrectBias (winner's-curse detection-bias correction) work

## Scope correction

The originally-discussed "model per-sample calibration in the likelihood" idea (a general
purity-miscalibration correction anchored on each sample's clonal peak) is **not** what
`data/CorrectBias.py` actually implements. Reading the real module clarified this: it implements
a **winner's-curse / detection-limit correction** specifically for low-CCF clusters --
`simulate_mutations` + `pre_compute_WCC` simulate the mutation-detection process at a sample's
real purity/coverage/local-CN profile to estimate how much a cluster's *apparent* CCF is biased
upward because its lowest-VAF mutations failed to get called at all (classic winner's-curse /
selection-bias truncation). This plan is scoped to what the module actually does. The broader
whole-sample, clonal-peak-anchored calibration idea is a separate, larger project (would need a
new joint per-sample nuisance parameter in the DP's likelihood, not just fixing this module) and
is not addressed here.

## Problem

`data/CorrectBias.py` (the "WC_Correction" module) was non-functional:
- **Blocked entirely under Python 3**: `muts_gen.next()` (Python 2 iterator protocol), bare
  `print` statements, and `[...] + range(...)` (Python 2 allowed list + range concatenation;
  Python 3's `range` is not a list) would all raise errors the moment the module was exercised.
- **A no-op even if it ran**: `simulate_mutations` skips any mutation with NaN local copy number
  (`if sum(np.isnan([a1, a2])) > 0: continue`), and `data/SomaticEvents.py`'s `clean_local_cn`
  had `float(cn1) if not np.nan else np.nan` -- since `np.nan` is truthy, `not np.nan` is always
  `False`, so this **unconditionally** set `local_cn_a1`/`local_cn_a2` to NaN regardless of the
  real value passed in. Every mutation would hit the NaN skip, so `simulate_mutations` would
  never simulate anything.
- **Unbounded convergence loop**: `pre_compute_WCC`'s inner `while abs(err) > 0.005 and
  cluster_pos > 0:` loop had no iteration cap -- if the simulated detection process oscillated
  instead of converging (plausible near `cluster_pos ~ 0`, where detection is sparse/discrete),
  this could run indefinitely.
- **Not wired up anywhere**: `pre_compute_WCC` was dead code, called from nowhere in the pipeline.

## Changes made

1. **`data/SomaticEvents.py`** (`clean_local_cn`): fixed the always-`False` conditional to
   actually check `np.isnan()` on the real (float-converted) input, so `local_cn_a1`/`local_cn_a2`
   are populated correctly instead of always being NaN. This also fixes `data/CorrectBias.py`'s
   consumption of these fields (previously always skipped) and, as a side effect, unblocks any
   other code that reads `mut.local_cn_a1`/`local_cn_a2` (e.g. output columns in
   `output/PhylogicOutput.py`).

2. **`data/CorrectBias.py`**:
   - `muts_gen.next()` -> `next(muts_gen)`.
   - bare `print` statements -> `logging.info`/`logging.warning`.
   - `[...] + range(...)` -> `[...] + list(range(...))` (a second, previously-unnoticed Python 3
     incompatibility found while porting this file -- `range` is a lazy view object in Python 3,
     not concatenable with a list).
   - `pre_compute_WCC` now takes a `max_iter` parameter (default 40) on its per-candidate-position
     convergence loop; if the cap is hit, the best estimate so far is kept and a warning is logged,
     instead of looping unboundedly.
   - New function `apply_wcc_correction(patient_data, cluster_ccfs, low_ccf_threshold=0.15,
     max_iter=40)`: fits a `pre_compute_WCC` correction curve per sample, then for every cluster
     whose per-sample CCF mode is below `low_ccf_threshold`, shifts that histogram's mass by the
     estimated correction delta (a rigid bin-shift + renormalize -- the reference algorithm only
     produces a point-position mapping, not a full re-derived posterior, so this is a
     transparent, matching-fidelity approximation) and renormalizes. Clusters/samples at or above
     the threshold are copied through unchanged. Only CCF *values* are adjusted; mutation-to-
     cluster assignments are untouched, matching the paper's documented scope.

3. **`output/PhylogicOutput.py`** (`write_patient_cluster_ccfs`): added an optional `filename`
   parameter (defaults to the existing `{indiv_name}.cluster_ccfs.txt` pattern, so existing
   callers are unaffected) so the corrected cluster CCFs can be written to a second file without
   duplicating the writer logic.

4. **`PhylogicNDT.py`**: two new CLI flags on the `Cluster` subparser:
   - `--correct_bias` (flag, off by default): run the correction after DP clustering.
   - `--correct_bias_ccf_threshold` (float, default `0.15`): the "low CCF" cutoff used to decide
     which cluster/sample entries get corrected.

5. **`Cluster/Cluster.py`**: after the normal cluster/mutation CCF output is written, if
   `--correct_bias` is set, calls `CorrectBias.apply_wcc_correction(...)` and writes the result to
   `{indiv_id}.cluster_ccfs.corrected.txt` via the now-parameterized `write_patient_cluster_ccfs`.
   This is purely additive -- normal output is unaffected whether or not the flag is used.

## What this does *not* do (explicitly out of scope for this branch)

- Does not touch mutation-to-cluster assignments -- only reported CCF values for already-finalized
  clusters.
- Does not feed corrected values into BuildTree automatically. A `--use_corrected_ccf` flag on
  BuildTree (reading `{indiv_id}.cluster_ccfs.corrected.txt` instead of the raw file) would be a
  natural, small follow-up once the correction itself has been validated.
- Does not parallelize the per-sample fitting loop (see Runtime section) -- each sample's
  `pre_compute_WCC` call runs sequentially in `apply_wcc_correction`.
- Does not address whole-sample purity-miscalibration/clonal-peak-anchored correction (see Scope
  correction above) -- only low-CCF, detection-limit-driven bias.

## Runtime impact as complexity scales up

- **Decoupled from the DP entirely.** `apply_wcc_correction` runs once, after clustering
  completes -- it has no interaction with `--n_iter`, mutation burden, or realized cluster count
  from the DP's perspective. Applying a fitted correction curve to a cluster's CCF histogram is a
  cheap interpolation + array shift, negligible per cluster.
- **Scales linearly with sample count, not with mutation burden.** Each sample gets one
  independent `pre_compute_WCC` fit; the simulation target inside it is a fixed 500 mutations
  regardless of how many real mutations your dataset has, so per-sample cost does not grow with
  mutation burden. But 10+ samples means 10+ of these fits, run sequentially in this
  implementation, so the total added wall-clock time is directly proportional to sample count.
- **Each individual fit is nontrivially expensive**: `pre_compute_WCC` evaluates 25 candidate CCF
  positions, each requiring a Newton-style convergence loop (now capped at `max_iter=40`), each
  iteration running 5 replicate 500-mutation simulations, each mutation requiring a 501-bin
  binomial pmf evaluation. This is real, bounded-but-nonzero added cost per sample.
- **Natural mitigation, not yet implemented**: since every sample's fit is fully independent, this
  loop is an ideal candidate for `concurrent.futures.ProcessPoolExecutor`, which would turn the
  linear-in-samples wall-clock cost into roughly "cost of the single slowest sample's fit" given
  enough cores. This wasn't implemented in this pass because `TumorSample`'s picklability across
  process boundaries (it holds hashtables, interval trees, and other nested state) hasn't been
  verified -- worth confirming before parallelizing, to avoid a fragile implementation that only
  sometimes works.
- **The `max_iter` safety cap directly bounds the worst case.** Before this fix, a single
  non-converging sample could make the whole correction step (and therefore the whole `Cluster`
  run, since it happens synchronously at the end) hang indefinitely. Now the worst case per sample
  per candidate position is bounded at 40 iterations x 5 simulations x 500-mutation draws --
  slower in the failure case, but never unbounded.

## Testing notes

Not yet validated against real or simulated data. Before relying on this:
- Confirm `clean_local_cn`'s fix doesn't change behavior anywhere it shouldn't (grep for
  `local_cn_a1`/`local_cn_a2` consumers -- e.g. `output/PhylogicOutput.py` columns -- and check
  the previously-always-NaN columns now populate sensibly).
- Run `--correct_bias` on a dataset with a known/expected low-CCF cluster and confirm the
  corrected CCF moves in the expected direction (down, since detection bias only inflates
  apparent CCF) and by a plausible magnitude.
- Verify `pre_compute_WCC`'s `max_iter` cap doesn't trigger routinely on typical data (if it does,
  40 may be too low, or there's a real convergence issue worth investigating rather than papering
  over with the cap).
