function [voxels_noneff, R2s_noneff, pv_noneff] = create_noneffcluster(X, dist_settings, hyparams)
    R2s_noneff = nan(1, hyparams.nb_noneffvx);
    pv_noneff = nan(1, hyparams.nb_noneffvx);

    voxels_noneff = dist_settings.vox_std_range .* randn(hyparams.nb_timesteps, ...
            hyparams.nb_noneffvx);
    for ith_voxel = 1: hyparams.nb_noneffvx
        [~, ~, ~, ~, STATS] = regress(voxels_noneff( :, ith_voxel), X);
        R2s_noneff(ith_voxel) = STATS(1);
        pv_noneff(ith_voxel) = STATS(3);
    end
end