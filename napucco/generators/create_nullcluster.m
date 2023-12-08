function [voxels_null, R2s_vxnull, pv_param_vxnull] = create_nullcluster(X, dist_settings, hyparams)
    R2s_vxnull = nan(1, hyparams.second_cluster);
    pv_param_vxnull = nan(1, hyparams.second_cluster);

    voxels_null = dist_settings.vox_std_range .* randn(hyparams.nb_timesteps, ...
            hyparams.second_cluster);
    for ith_voxel = 1: hyparams.second_cluster
        [~, ~, ~, ~, STATS] = regress(voxels_null( :, ith_voxel), X);
        R2s_vxnull(ith_voxel) = STATS(1);
        pv_param_vxnull(ith_voxel) = STATS(3);
    end
end