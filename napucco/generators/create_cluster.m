function [voxels, betas, R2s, pv_param] = create_cluster(X, R2_target, nb_voxels, logger, opt)
%CREATE_CLUSTER creates a cluster of N voxels each with regards to one
%subject rapresented as a matrix which columns are features and the rows the
%observations.
% [voxels, betas, R2s] = create_cluster(X, R2distro, nb_voxels)  with respect to
% the feature matrix X, the function returns a voxel matrix with 
% nb_voxels colums corresponding to each voxel and size(X,1) rows
% corresponding to each observation. The desired R2 can be given as a pair of mean
% and std deviation.
    arguments
        X (:, :) double
        R2_target (1,:) double
        nb_voxels int32
        logger
        opt.R2distro (1,:) double
        opt.R2_tollerance double = 0.01
        opt.rnd_is_normal int32 = 1
        opt.vx_std_range = 0.1 %stdev of noise applied to voxel activity
        opt.vox_variability = 0 %within-subject sampling from same dist of R2
    end

    voxels =  nan(size(X, 1), nb_voxels);
    betas =  nan(size(X, 2), nb_voxels);
    R2s =  nan(1, nb_voxels);
    pv_param =  nan(1, nb_voxels);

    parfor i = 1:nb_voxels
        R2_target_final = R2_target;
        if opt.vox_variability > 0
            R2_target_final = 0;
            while R2_target_final<=0.005 || R2_target_final>0.8
                R2_target_final = R2_target - abs(opt.vox_variability*randn(1));
            end
        end
        
        [Y, ~, R2, ~] = createSingVoxBehav( X, R2_target_final, ...
            opt.R2_tollerance, opt.vx_std_range, logger);
        voxels( :, i) = Y;
        R2s(i) = R2;
    end
    if opt.vox_variability > 0 % if there is variance within subject (fpr clusterbased)
        sorted_indices = cubic_radial_sort(nb_voxels);
        [R2s, outdata_indices] = sort(R2s, "descend");
        R2s = R2s(sorted_indices); 
        voxels = voxels(:, outdata_indices(sorted_indices));
    end
end

