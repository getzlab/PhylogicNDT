##WC_Correction

import logging
import random
import itertools
import numpy as np
from scipy.stats import binom


def est_ccf_dist(alt, ref, allele_cn, PURITY, other_cn, mult):
    ccf_dist = np.zeros(501)
    ccf_space = np.linspace(0, 1, 501)

    af_points = []

    for mult_1_bin_val_idx, mult_1_bin_val in enumerate(ccf_space):
        af = (mult * PURITY * mult_1_bin_val) / (
                (allele_cn * PURITY + other_cn * PURITY + 2 * (1 - PURITY)) + (
                PURITY * mult + (allele_cn - mult) * PURITY - allele_cn * PURITY) * (mult_1_bin_val))
        af_points.append(af)
        # ccf_dist[mult_1_bin_val_idx] = m1_draw

    ccf_dist[:] = binom._pmf(alt, alt + ref, np.array(af_points))

    ccf_dist[np.isnan(ccf_dist)] = 0.
    return ccf_dist / sum(ccf_dist)


def simulate_mutations(PURITY, NMUT, ccf, mutations):
    random.shuffle(mutations)
    muts_gen = itertools.cycle(mutations)

    clust_dist = np.zeros(501)

    mut_num = 0
    total_mut = 0
    while mut_num < NMUT and total_mut < 100 * NMUT:

        mut = next(muts_gen)
        if mut.type != "SNP":
            continue
        real_cov = mut.alt_cnt + mut.ref_cnt
        if real_cov == 0:
            continue
        a1, a2 = mut.local_cn_a1, mut.local_cn_a2
        if sum(np.isnan([a1, a2])) > 0:
            continue
        cn = a1 + a2

        if cn == 0:
            continue

        mult = np.random.choice([x for x in [a1, a2] if x > 0] + list(range(1, int(max([a1, a2])) + 1)))

        af = (ccf * mult * PURITY) / (
                float(ccf) * mult * PURITY + ccf * (cn - mult) * PURITY + (1 - ccf) * (cn) * PURITY + 2 * (
                1.0 - PURITY))

        alt_count = np.random.binomial(real_cov, af)
        ref_count = int(real_cov - alt_count)

        total_mut += 1

        if alt_count >= 1:
            mut_num += 1
            ccf_dist = est_ccf_dist(alt_count, ref_count, a1, PURITY, a2, mult)
            clust_dist += np.log(ccf_dist)

    return clust_dist


def pre_compute_WCC(sample, max_iter=40):
    """
    Fit a winner's-curse / detection-bias correction curve for one sample: for a grid of candidate
    "observed" cluster CCF positions, iteratively simulate the mutation-detection process at this
    sample's real purity/coverage/local-CN profile (via simulate_mutations) to find the "true"
    cluster CCF that would, after detection bias from low-VAF mutations failing to be called,
    appear at that observed position. Returns an interp1d mapping observed -> corrected CCF.

    max_iter caps the per-candidate-position convergence loop, which previously had no bound and
    could run indefinitely if the simulated detection process oscillated instead of converging
    (most likely near cluster_pos ~ 0, where detection is sparse/discrete). If the cap is hit, the
    best estimate found so far is kept and a warning is logged, rather than a silent/unbounded hang.
    """
    logging.info("Pre-computing winner's-curse correction (WCC) for sample {}".format(
        getattr(sample, 'sample_name', '?')))
    corr_y = []

    for obs_cluster_pos in list(np.logspace(0, 1, 25) / 10. - 0.1) + [1.]:

        err = 1
        cluster_pos = obs_cluster_pos
        n_iter = 0
        while abs(err) > 0.005 and cluster_pos > 0:
            if n_iter >= max_iter:
                logging.warning(
                    "WCC fit for sample {} did not converge within {} iterations at observed "
                    "position {:.3f} (last err={:.4f}); using best estimate so far.".format(
                        getattr(sample, 'sample_name', '?'), max_iter, obs_cluster_pos, err))
                break
            sim_res = np.average(
                [np.argmax(simulate_mutations(sample.purity, 500, cluster_pos, sample.concordant_variants)) / 500. for x
                 in range(5)])
            err = obs_cluster_pos - sim_res  # distance between observed cluster position, and the observed simulated cluster position.
            cluster_pos += err / 2.
            n_iter += 1
        corr_y.append(max(cluster_pos, 0))

    from scipy.interpolate import interp1d

    logging.info("Done pre-computing WCC for sample {}".format(getattr(sample, 'sample_name', '?')))

    return interp1d(list(np.logspace(0, 1, 25) / 10. - 0.1) + [1.], corr_y)


def apply_wcc_correction(patient_data, cluster_ccfs, low_ccf_threshold=0.15, max_iter=40):
    """
    Apply the winner's-curse/detection-bias correction to cluster CCF histograms whose per-sample
    mode falls below low_ccf_threshold. Mutations that fail to be called at very low allele
    fractions bias a low-CCF cluster's apparent CCF upward; this fits a per-sample correction curve
    (pre_compute_WCC) from that sample's real purity/coverage/local-CN profile and shifts affected
    histograms back down toward their estimated true position.

    Only cluster CCF *values* are adjusted -- mutation-to-cluster assignments are untouched, matching
    PhylogicNDT CorrectBias's documented scope (paper: "adjustment of clustered CCF values and
    cluster sizes... in subclones with low CCF").

    Args:
        patient_data: data.Patient instance (already clustered; sample_list's concordant_variants
            must be populated, i.e. called after preprocess_samples()/run_DP_ND()).
        cluster_ccfs: dict mapping cluster_id -> list of per-sample CCF histograms (as built in
            Cluster/Cluster.py from ClusteringResults.clust_CCF_dens).
        low_ccf_threshold: only histograms whose mode is below this CCF are corrected; clusters/
            samples at or above it are copied through unchanged.
        max_iter: passed to pre_compute_WCC's per-position convergence loop safety cap.

    Returns:
        A new dict, same shape as cluster_ccfs, with corrected histograms substituted in for
        low-CCF entries.
    """
    grid_size = patient_data.sample_list[0].ccf_grid_size

    # Each sample's correction curve is independent of every other sample's -- this loop is an
    # embarrassingly parallel candidate (e.g. concurrent.futures.ProcessPoolExecutor mapping
    # pre_compute_WCC over patient_data.sample_list) if this becomes a bottleneck with many samples.
    # Left sequential here since TumorSample's picklability across process boundaries hasn't been
    # verified; see PLAN_correctbias.md.
    correction_fns = {}
    for sample in patient_data.sample_list:
        correction_fns[sample.sample_name] = pre_compute_WCC(sample, max_iter=max_iter)

    corrected = {}
    for cluster_id, per_sample_hists in cluster_ccfs.items():
        corrected_hists = []
        for sample_index, sample in enumerate(patient_data.sample_list):
            hist = np.array(per_sample_hists[sample_index], dtype=np.float64)
            observed_pos = np.argmax(hist) / float(grid_size - 1)

            if observed_pos >= low_ccf_threshold:
                corrected_hists.append(hist)
                continue

            fn = correction_fns[sample.sample_name]
            try:
                corrected_pos = float(fn(observed_pos))
            except ValueError:
                # observed_pos fell outside the fitted interpolation domain; leave uncorrected
                # rather than extrapolate.
                corrected_hists.append(hist)
                continue

            shift_bins = int(round((corrected_pos - observed_pos) * (grid_size - 1)))
            if shift_bins == 0:
                corrected_hists.append(hist)
                continue

            # Shift the histogram's mass by the estimated correction delta and renormalize. This is
            # a rigid-shift approximation of the correction (the reference algorithm only produces a
            # point-position mapping, not a full re-derived posterior), applied consistently in both
            # directions.
            shifted = np.roll(hist, shift_bins)
            if shift_bins > 0:
                shifted[:shift_bins] = 0.
            else:
                shifted[shift_bins:] = 0.
            total = shifted.sum()
            if total > 0:
                shifted = shifted / total
            else:
                shifted = hist  # degenerate shift; fall back to uncorrected rather than an all-zero row
            corrected_hists.append(shifted)

        corrected[cluster_id] = corrected_hists

    return corrected
