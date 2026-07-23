# Plan: Pseudo-count reparametrization of the DP's cluster-count prior (`--Pi_k_strength`)

## Problem

`--Pi_k_mu`/`--Pi_k_r` (the negative-binomial prior mean/dispersion over the number of DP
clusters, K) often have little visible effect on the realized clustering, especially with many
(10+) samples and high mutation burden.

This is intentional in one sense: per the lab's own design notes, alpha (the DP concentration
parameter) is deliberately re-*learned* from the data every iteration rather than held fixed, so
that the resulting clustering reflects the actual underlying data rather than an arbitrarily
imposed prior belief about K. That part of the design is sound and shouldn't change.

The problem is narrower: `Pi_k_mu`/`Pi_k_r` only shape the **mean and variance** of the target
K-distribution used to derive the alpha prior's `(a, b)` Gamma parameters
(`Cluster/DpEngine.py`, `get_gamma_prior_from_k_prior`/`init_dp_prior`). Nothing controls how
*strongly* that belief is held, independent of its shape. Each iteration, alpha is resampled via
the Escobar & West (1995) conjugate update (`sample_gamma_cond_N_k`), where (informally):

```
shape ~ a + realized_k
rate  ~ b + log(N)
```

`a` and `b` are small, fixed constants from the initial fit; `realized_k` (the chain's current
cluster count) and `log(N)` (mutation burden) are live, data-driven quantities of comparable or
larger magnitude once N/k grow. So the fixed `(a, b)` becomes a proportionally smaller share of
the resampled alpha as mutation burden and realized cluster count grow -- which is exactly why
tuning `Pi_k_mu`/`Pi_k_r` has little visible effect on complex, high-sample-count datasets: the
belief they encode gets diluted by data volume, even though nothing about the user's stated
belief changed.

## Fix

Add a **pseudo-count / prior-strength** multiplier that scales the *already-fitted* `(a, b)` up
by a constant factor, without touching the fitting logic or the resampling math at all.

For a Gamma(a, b) distribution: mean = a/b, variance = a/b^2. Scaling both `a` and `b` by the
same constant `c` leaves the mean **exactly unchanged** (so the belief about K's central value is
untouched) while shrinking the variance by `1/c` -- precisely the standard pseudo-count
interpretation: "the same central belief, now backed by `c` times as much confidence." This keeps
alpha fully resampled and data-adaptive every iteration (nothing is frozen) -- it's just harder
for the realized `k`/`log(N)` terms to drag it away from the stated prior.

## Changes made

1. **`PhylogicNDT.py`**: new CLI flag `--Pi_k_strength` on the `Cluster` subparser, float,
   default `1.0` (reproduces current behavior exactly).
2. **`Cluster/Cluster.py`**: passes `args.Pi_k_strength` into `ClusterEngine.run_DP_ND(...)`.
3. **`Cluster/ClusterEngine.py`**: `run_DP_ND` accepts `prior_strength` and passes it to the
   `DpEngine` constructor. (Also fixed the pre-existing dead default `N_iter=5` -- vs. the CLI's
   real default of 250, always passed explicitly -- to `N_iter=250`, removing a landmine for any
   future direct caller of this method.)
4. **`Cluster/DpEngine.py`**: `DpEngine.__init__` accepts `prior_strength` (default `1.0`) and
   applies it immediately after the existing fit:
   ```python
   DP_prior = init_dp_prior(len(data._hist_array), Pi_k)
   Pi_gamma_a = DP_prior["a"] * prior_strength
   Pi_gamma_b = DP_prior["b"] * prior_strength
   ```
   `get_gamma_prior_from_k_prior`, `init_dp_prior`, and `sample_gamma_cond_N_k` are all
   **unmodified** -- this is a pure post-hoc rescale of their output.

## How to use it

Start at the default (`1.0`, no change). On a dataset where `--Pi_k_mu`/`--Pi_k_r` currently seem
to have no effect, increase `--Pi_k_strength` (try 5, 20, 100...) until the realized number of
clusters visibly responds. Because the mean is preserved exactly by construction, this is a much
more predictable knob to turn than re-guessing `Pi_k_mu` itself -- you're only ever adjusting how
firmly the belief is held, never what it's centered on.

## Runtime impact as complexity scales up

- **Direct cost: effectively zero.** The change is a single multiplication applied once, at
  `DpEngine` initialization, before the MCMC loop starts. It does not change the number of Gibbs
  iterations, the per-iteration cost (still driven by `n_mutations x n_clusters x n_samples x
  num_bins`), or the cost of the `(a, b)` fitting procedure itself (the 5x5-grid Nelder-Mead
  search already runs today regardless of `prior_strength`).
- **Indirect effect: likely a net *reduction* in downstream cost for high-complexity datasets.**
  By keeping alpha anchored nearer the stated `Pi_k_mu` even as `k`/`N` grow, a well-tuned
  `--Pi_k_strength` should suppress some of the spurious cluster-count inflation that otherwise
  compounds in the 10+-sample / high-mutation-burden regime. Fewer realized clusters means:
  - less per-DP-iteration cost (which scales with cluster count),
  - a better chance the K-selection stability check (the "at least 10% of iterations at this K"
    rule) converges within the existing `--n_iter` budget instead of landing on a diffuse
    posterior over K,
  - fewer clusters handed to BuildTree, whose per-sweep move space grows roughly quadratically
    (or worse) with cluster count.
  So, unlike a change that adds new computation, this one plausibly makes complex runs *faster*
  end-to-end, not just no-slower.
- **Caveat (accuracy, not runtime):** if `--Pi_k_strength` is set too high and the data genuinely
  supports more clusters than `Pi_k_mu` implies, real splits could be suppressed along with
  spurious ones. This is a tuning/validation question (calibrate on a dataset with a known/
  expected cluster count), not a performance concern.

## Testing notes

Not yet validated against real or simulated data. Before relying on this:
- Confirm on a known high-N/high-sample dataset that increasing `--Pi_k_strength` measurably
  shifts realized K back toward `Pi_k_mu`, where today changing `Pi_k_mu`/`Pi_k_r` alone does not.
- Sanity-check that `--Pi_k_strength=1.0` reproduces bit-identical (or statistically
  indistinguishable, given the MCMC is stochastic) results to the current `master` behavior on at
  least one existing test case, given the same `--seed`.
