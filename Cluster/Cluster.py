import os
import logging


# Main run method

def run_tool(args):
    print(args)
    if not args.maf and not args.html:
        print("No output type specified, specify any combination of --html and --maf")

    # sys.path.append(args.phylogic_path)
    # data storage objects
    import data.Patient as Patient
    import data.Sample as Sample

    import ClusterEngine
    import DpEngine
    if not os.path.isfile(args.PoN):
        logging.warning("PanelofNormals (PoN) inaccessible or not specified - not using PoN!")
        PoN = False
    else:
        PoN = args.PoN

    # init a Patient
    patient_data = Patient.Patient(artifact_blacklist=args.artifact_blacklist,
                                   PoN_file=PoN, indiv_name=args.indiv_id, artifact_whitelist=args.artifact_whitelist,
                                   min_coverage=args.min_cov, use_indels=args.use_indels,
                                   impute_missing=args.impute_missing,
                                   driver_genes_file=args.driver_genes_file)

    # delete_auto_bl=args.Delete_Blacklist,
    # Load sample data

    if args.sif:  # if sif file is specified
        sif_file = open(args.sif)

        for file_idx, file_line in enumerate(sif_file):

            ##for now, assume input file order is of the type sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint
            if not file_idx:
                continue
            if file_line.strip('\n').strip == "":
                continue  # empty rows
            smpl_spec = file_line.strip('\n').split('\t')
            sample_id = smpl_spec[0]
            maf_fn = smpl_spec[1]
            seg_fn = smpl_spec[2]
            purity = float(smpl_spec[3])
            timepoint = float(smpl_spec[4])
            print(timepoint)
            patient_data.addSample(maf_fn, sample_id, timepoint_value=timepoint, grid_size=args.grid_size,
                                   _additional_muts=None, seg_file=seg_fn,
                                   purity=purity, input_type=args.maf_input_type)

            if len(patient_data.sample_list) == args.n_samples:  # use only first N samples
                break

    else:  # if sample names/files are specified directly on cmdline

        # sort order on timepoint or order of entry on cmdline if not present
        print(args.sample_data)
        for idx, sample_entry in enumerate(args.sample_data):
            ##for now, assume input order is of the type sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint
            smpl_spec = sample_entry.strip('\n').split(':')
            sample_id = smpl_spec[0]
            maf_fn = smpl_spec[1]
            seg_fn = smpl_spec[2]
            purity = float(smpl_spec[3])
            timepoint = float(smpl_spec[4])

            patient_data.addSample(maf_fn, sample_id, timepoint_value=timepoint, grid_size=args.grid_size,
                                   _additional_muts=None,
                                   seg_file=seg_fn,
                                   purity=purity)
            if len(patient_data.sample_list) == args.n_samples:  # use only first N samples
                break

    patient_data.get_arm_level_cn_events()
    patient_data.preprocess_samples()

    if getattr(args, 'remove_private_muts', False):
        patient_data.flag_private_mutations(present_cutoff=args.private_ccf_cutoff,
                                            absent_cutoff=args.private_absent_cutoff)

    # TODO: how 1D (one sample) is handeled
    DP_Cluster = ClusterEngine.ClusterEngine(patient_data)  # html_out=args.indiv_id + '.html')

    DP_Cluster.run_DP_ND(N_iter=args.iter, PriorK={'r': args.Pi_k_r, 'mu': args.Pi_k_mu},
                         mode="maxpear", seed=args.seed)

    patient_data.ClusteringResults = DP_Cluster.results
    patient_data.cluster_temp_removed()
    # sample_id = args.indiv_id

    # Output and visualization
    import output.PhylogicOutput
    phylogicoutput = output.PhylogicOutput.PhylogicOutput()
    cluster_ccfs = {i + 1: ccf_hists for i, ccf_hists in enumerate(patient_data.ClusteringResults.clust_CCF_dens)}
    phylogicoutput.write_patient_cluster_ccfs(patient_data, cluster_ccfs)
    phylogicoutput.write_patient_mut_ccfs(patient_data, cluster_ccfs)
    phylogicoutput.write_patient_cnvs(patient_data, cluster_ccfs)
    phylogicoutput.write_patient_unclustered_events(patient_data)
    phylogicoutput.plot_1d_clusters('{}.cluster_ccfs.txt'.format(patient_data.indiv_name))
    phylogicoutput.plot_1d_mutations('{}.mut_ccfs.txt'.format(patient_data.indiv_name))

    if getattr(args, 'correct_bias', False):
        import data.CorrectBias as CorrectBias
        corrected_cluster_ccfs = CorrectBias.apply_wcc_correction(
            patient_data, cluster_ccfs, low_ccf_threshold=args.correct_bias_ccf_threshold)
        phylogicoutput.write_patient_cluster_ccfs(
            patient_data, corrected_cluster_ccfs,
            filename='{}.cluster_ccfs.corrected.txt'.format(patient_data.indiv_name))

    if not args.buildtree:  # run only Clustering tool
        phylogicoutput.generate_html_from_clustering_results(patient_data.ClusteringResults, patient_data,
                                                             drivers=patient_data.driver_genes,
                                                             treatment_file=args.treatment_data)

    else:  # run build tree next
        import BuildTree.BuildTree
        args.cluster_ccf_file = '{}.cluster_ccfs.txt'.format(patient_data.indiv_name)
        args.mutation_ccf_file = '{}.mut_ccfs.txt'.format(patient_data.indiv_name)
        args.n_iter = args.iter
        # Auto-blacklist any dedicated private-mutation clusters (see flag_private_mutations /
        # _assign_private_clusters) from BuildTree, in addition to whatever the user passed in.
        # NOTE: this used to unconditionally reset to None, discarding any user-supplied
        # --blacklist_cluster ids when chaining Cluster -> BuildTree in one invocation.
        # NOTE 2: --blacklist_cluster is only defined as a CLI arg on the BuildTree subparser, not
        # on Cluster's -- so when running the combined Cluster+BuildTree flow (-rb), args never
        # has a blacklist_cluster attribute at all (not just None), and a direct args.blacklist_cluster
        # read raised AttributeError. Use getattr with a default instead.
        private_cluster_ids = [str(c) for c in getattr(DP_Cluster.results, 'private_cluster_ids', [])]
        args.blacklist_cluster = (getattr(args, 'blacklist_cluster', None) or []) + private_cluster_ids
        BuildTree.BuildTree.run_tool(args)
