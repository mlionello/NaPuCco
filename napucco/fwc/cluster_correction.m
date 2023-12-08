function results = cluster_correction(voxels, subj_range, pies_thr, opts, logger, cdt, parametric, non_eff_range, padding_range)
arguments
    voxels;
    subj_range;
    pies_thr;
    opts;
    logger;
    cdt;
    parametric = 1;
    non_eff_range = [];
    padding_range = [];
end
    results = struct();
    voxels_eff = permute(voxels.eff, [2 1 3]); % to match [nb_perm + 1, nb_subj, nb_voxels]
    voxels_null = permute(voxels.null, [2 1 3]);
    
    nb_perm = size(voxels_eff, 1);
    nb_rangesSubj = size(voxels_eff, 2);
    nb_effv = size( voxels_eff, 3);
    nb_noneffv = size( voxels_null, 3);
    
    % cluster formation threshold is applied
    if parametric
        voxels_eff_pv = 1 - chi2cdf(voxels_eff, repmat(double(2*subj_range), ...
            size(voxels_eff, 1), 1, size(voxels_eff, 3)));
        
        voxels_null_pv = 1 - chi2cdf(voxels_null, repmat(double(2*subj_range), ...
            size(voxels_null, 1), 1, size(voxels_null, 3))); %subj_range is 1 by N
    else
        voxels_eff_pv = tiedrank(-voxels_eff) / size(voxels_eff,1); % MINUS IF FISHER
        voxels_null_pv = tiedrank(-voxels_null) / size(voxels_null,1); % MINUS IF FISHER
    end
    cft_eff = voxels_eff_pv < cdt;
    cft_noneff = voxels_null_pv < cdt;

    if opts.power
        [K_power, numb_effcftvox, K_effcftvox, m, n, padding_range, power] = clust_corr_power( ...
            cft_eff, cft_noneff, nb_effv, nb_noneffv, nb_perm, ...
            nb_rangesSubj,pies_thr, padding_range, logger);
        results.K_power = K_power;
        results.numb_effcftvox = numb_effcftvox;
        results.K_effcftvox = K_effcftvox;
        results.m = m;
        results.n = n;
        results.padding_range = padding_range;
        results.power = power;
    end
    if opts.fpr
        [fpr, non_eff_range, K_fpr] = clus_corr_fpr(cft_noneff, nb_noneffv, nb_rangesSubj, nb_perm, pies_thr, non_eff_range, logger);
        results.fpr = fpr;
        results.non_eff_range = non_eff_range;
        results.K_fpr = K_fpr;
    end
    
    results.pies_thr = pies_thr;
    results.nb_perm = nb_perm ;
    results.subj_range = subj_range;
    results.nb_effv = nb_effv;
    results.nb_noneffv = nb_noneffv;
    results.cdt = cdt;
end

