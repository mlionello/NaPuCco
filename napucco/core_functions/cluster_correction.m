function results = cluster_correction(voxels, subj_range, pies_thr, opts, logger, cdt, len_pnoneff, non_eff_range, padding_range, parametric)
arguments
    voxels;
    subj_range;
    pies_thr;
    opts;
    logger;
    cdt;
    len_pnoneff;
    non_eff_range = [];
    padding_range = [];
    parametric = 1;
end
    results = struct();
    eff = permute(voxels.eff, [2 1 3]); % to match [nb_perm + 1, nb_subj, nb_voxels]
    noneff = permute(voxels.noneff, [2 1 3]);
    
    nb_perm = size(eff, 1);
    nb_rangesSubj = size(eff, 2);
    nb_effv = size( eff, 3);
    nb_noneffv = size( noneff, 3);
    
    % cluster formation threshold is applied
    logger.println("calculating chi2df...")
    [eff_pv, noneff_pv] = get_chi2_dist_from_cluster(eff, ...
        noneff, parametric, subj_range, len_pnoneff);
    cft_eff = eff_pv < cdt;
    cft_noneff = noneff_pv < cdt;

    % calculating experiments passing thresholds for power calculation:
    if opts.power
        logger.println("calculating power...")
        [K_max, numb_effcftvox, K_effcftvox, m, n, padding_range, passing_esperiment] = clust_corr_fn( ...
            cft_eff, cft_noneff, nb_effv, nb_noneffv, nb_perm, ...
            nb_rangesSubj, pies_thr, padding_range, logger);
        results.K_tfp_max = K_max;
        results.numb_effcftvox = numb_effcftvox;
        results.K_effcftvox = K_effcftvox;
        results.m = m;
        results.n = n;
        results.padding_range = padding_range;
        results.tfp = passing_esperiment;
    end
    
    % calculating experiments passing thresholds for fpr calculation:
    if opts.fpr
        logger.println("calculating fpr...")
        [passing_esperiment, non_eff_range, K_fpr] = clus_corr_fp(cft_noneff, nb_noneffv, nb_rangesSubj, nb_perm, pies_thr, non_eff_range);
        results.fn = passing_esperiment;
        results.non_eff_range = non_eff_range;
        results.K_fn_max = K_fpr;
    end
    
    results.pies_thr = pies_thr;
    results.nb_perm = nb_perm ;
    results.subj_range = subj_range;
    results.nb_effv = nb_effv;
    results.nb_noneffv = nb_noneffv;
    results.cdt = cdt;
end

