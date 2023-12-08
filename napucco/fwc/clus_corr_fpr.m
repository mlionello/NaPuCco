function [fpr, non_eff_range, K_fpr] = clus_corr_fpr(cft_noneff, nb_noneffv, nb_rangesSubj, nb_perm, pies_thr, non_eff_range, logger)
max_non_eff_cubic_size = cubicrt(nb_noneffv);
if isempty(non_eff_range)
    non_eff_range = 1: max_non_eff_cubic_size;
else
    non_eff_range = non_eff_range(non_eff_range<max_non_eff_cubic_size);
end

K_fpr = nan( nb_perm, nb_rangesSubj, length(non_eff_range));

parfor j = 1: nb_perm
    %logger.print("; perm %d/%d\r", j, size( voxels_eff_pv, 1) );

    for n_subj = 1: nb_rangesSubj
        % cluster formation threshold is applied to the non-effected voxels
        cft_noneff_jth = squeeze(cft_noneff(j, n_subj, :));

        % FPR
        for null_width = non_eff_range
            k = get_k(cft_noneff_jth(1: null_width^3));
            K_fpr( j, n_subj, null_width) = max(k);
        end
    end
end
logger.clear();

K_size = size(K_fpr);
K_tied = tiedrank( -K_fpr) / nb_perm;
for p_i = 1: length(pies_thr)
    fpr(1: K_size(2), 1: K_size(3), p_i) = K_tied(1, :, :) < pies_thr(p_i);
end
non_eff_range = non_eff_range;
K_fpr = K_fpr;
end
