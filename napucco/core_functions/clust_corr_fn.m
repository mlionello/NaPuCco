function [K_max, numb_effcftvox, K_effcftvox, m, n, padding_range, true_positives] = clust_corr_fn( cft_eff, cft_noneff, nb_effv, nb_noneffv, nb_perm, nb_rangesSubj, pies_thr, padding_range, logger)
    [n, ~, eff_width, m] = resizecubicpadding( ...
        nb_effv, nb_noneffv, logger);
    if isempty(padding_range)
        padding_range =  1: m;
    else
        padding_range = padding_range(padding_range <= m);
    end
    nb_padding_range = numel(padding_range);
    
    K_max = nan( nb_perm, nb_rangesSubj, nb_padding_range);
    K_effcftvox = nan( nb_perm, nb_rangesSubj, nb_padding_range);
    numb_effcftvox = nan( nb_perm, nb_rangesSubj, nb_padding_range);
    
    tic;
    parfor_progress(nb_perm);
    parfor (j = 1: nb_perm, min(8, gcp("nocreate").NumWorkers))
        parfor_progress;
    
        for n_subj = 1: nb_rangesSubj
            % cluster formation threshold is applied to the non-effected voxels
            cft_noneff_jth = squeeze(cft_noneff(j, n_subj, :));
    
            cft_eff_jth = squeeze(cft_eff(j, n_subj, :));
            % applying different padding sizes of non effected voxels
            for w = 1 : nb_padding_range
                tot_noneff_vx = get_nb_noneff_vx(eff_width, padding_range(w));
                flattenpadding = encubator( ...
                    cft_eff_jth, ...
                    cft_noneff_jth(1: tot_noneff_vx), ...
                    eff_width, ...
                    padding_range(w));
    
                k = get_k(flattenpadding);
                K_max( j, n_subj, w) = max(k);
    
                numb_effcftvox(j, n_subj, w) = sum(cft_eff_jth);
                K_effcftvox(j, n_subj, w) = max(get_k(cft_eff_jth));
            end
        end
    end
    parfor_progress(0);
    toc;
    
    K_size = size(K_max);
    % comparsion between null distribution of permutation and
    % the target Ks from voxels_eff_ranges
    K_tied = tiedrank( -K_max) / nb_perm;
    % For the results from the null distro tied rank comparison, if the
    % the results rank below the p_threhsold used in the mask above, the
    % expriemtns passes, it failes otherwise.
    for p_i = 1: length(pies_thr)
        true_positives( 1 : size(K_max, 2), 1 : size(K_max, 3), p_i) = K_tied(1, :, :) < pies_thr(p_i);
    end

end
