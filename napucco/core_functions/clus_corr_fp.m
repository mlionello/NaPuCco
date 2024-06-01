function [false_positives, non_eff_range, K_max] = clus_corr_fp(cft_noneff, ...
    nb_noneffv, nb_rangesSubj, nb_perm, pies_thr, non_eff_range)

max_non_eff_cubic_size = cubicrt(nb_noneffv);
if isempty(non_eff_range)
    non_eff_range = 1 : max_non_eff_cubic_size;
else
    non_eff_range = non_eff_range(non_eff_range <= max_non_eff_cubic_size);
end

K_max = nan( nb_perm, nb_rangesSubj, length(non_eff_range));
nb_non_eff_range = numel(non_eff_range);
tic;
parfor_progress(nb_perm);
parfor (j = 1: nb_perm, 8)
    parfor_progress;

    for n_subj = 1: nb_rangesSubj
        cft_noneff_jth = squeeze(cft_noneff(j, n_subj, :));

        for w = 1:nb_non_eff_range
            K = get_k(cft_noneff_jth(1 : non_eff_range(w) ^ 3));
            K_max( j, n_subj, w) = max(K);
        end
    end
end
parfor_progress(0);
toc;

K_tied = tiedrank( -K_max) / nb_perm;
for p_i = 1: length(pies_thr)
    false_positives(:, :, p_i) = squeeze(K_tied(1, :, :) < pies_thr(p_i));
end
end
