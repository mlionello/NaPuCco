function [valid_indices] = resample_nullR2(R2snull, stdev_sample_fromnull, nb_noneffvoxels)
%RESAMPLE_NULLR2 
% R2snull := hyparams.nb_subj by hyparams.second_cluster
nb_subj = size(R2snull,1);
valid_indices = nan(nb_subj, nb_noneffvoxels);
for subj_i = 1:nb_subj
    init_mean = mean(R2snull(subj_i, :));
    init_std = std(R2snull(subj_i, :));
    targetR2 = init_mean + stdev_sample_fromnull*randn(nb_noneffvoxels, 1);
    for r2_i = 1: nb_noneffvoxels
        [~, valid_indices(subj_i, r2_i)] = min(abs(R2snull(subj_i, :) - targetR2(r2_i)));
    end
end
end
