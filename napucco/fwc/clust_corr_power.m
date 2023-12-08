function [K_power, numb_effcftvox, K_effcftvox, m, n, padding_range, power] = clust_corr_power( cft_eff, cft_noneff, nb_effv, nb_noneffv, nb_perm, nb_rangesSubj, pies_thr, padding_range, logger)
[n, padding_indices, eff_width, m] = resizecubicpadding( ...
    nb_effv, nb_noneffv, logger);
K_power = nan( nb_perm, nb_rangesSubj, m);
if isempty(padding_range)
    padding_range =  1: m;
else
    padding_range = padding_range(padding_range<=m);
end

for j = 1: nb_perm
    %logger.print("; perm %d/%d\r", j, size( voxels_eff_pv, 1) );

    for n_subj = 1: nb_rangesSubj
        % cluster formation threshold is applied to the non-effected voxels
        cft_noneff_jth = squeeze(cft_noneff(j, n_subj, :));

        % POWER
        cft_eff_jth = squeeze(cft_eff(j, n_subj, :));
        % applying different padding sizes of non effected voxels
        for null_width = padding_range
            tot_null_vx = get_nb_null_vx(eff_width, null_width);
            flattenpadding = encubator( ...
                cft_eff_jth, ...
                cft_noneff_jth(1: tot_null_vx), ...
                eff_width, ...
                null_width);

            k = get_k(flattenpadding);
            K_power( j, n_subj, null_width) = max(k);

            numb_effcftvox(j, n_subj, null_width) = sum(cft_eff_jth);
            K_effcftvox(j, n_subj, null_width) = max(get_k(cft_eff_jth));
        end
    end
end
logger.clear();
logger.flush();

K_size = size(K_power);
% comparsion between null distribution of permutation and
% the target Ks from voxels_eff_ranges
K_tied = tiedrank( -K_power) / nb_perm;
% For the results from the null distro tied rank comparison, if the
% the results rank below the p_threhsold used in the mask above, the
% expriemtns passes, it failes otherwise.
for p_i = 1: length(pies_thr)
    power( 1:K_size(2), 1:K_size(3), p_i) = K_tied(1, :, :) < pies_thr(p_i);
end

end

