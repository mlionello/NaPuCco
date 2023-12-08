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

%     if size(R2distro, 2) == 1
%         if ~ opt.rnd_is_normal
%             error('missing ending range value for uniform distribution of R2 in cluster')
%         end
%         %R2distro = [R2distro, 0.02]; %R2 variance of voxels activity across subjects
%     end
    voxels =  nan(size(X, 1), nb_voxels);
    betas =  nan(size(X, 2), nb_voxels);
    R2s =  nan(1, nb_voxels);
    pv_param =  nan(1, nb_voxels);

    %add here R2 sampling if i want smoothed voxels activities  in the clusters
%     if ~opt.rnd_is_normal
%         R2_target = R2distro(1) + 0.1.*randn(nb_voxels, 1);
%         R2_target = max(min(R2distro, 1), 0);
%     else
%         R2_target = unifrnd(R2distro(1), R2distro(2), nb_voxels, 1);
%     end
    parfor i = 1:nb_voxels
        %logger.print("; voxel: %d/%d\r", i, nb_voxels);
        R2_target_final = R2_target;
        if opt.vox_variability > 0
            R2_target_final = 0;
            while R2_target_final<=0.005 || R2_target_final>0.8
%                 a = ((1 - R2_target) / opt.vox_variability^2 - 1 / R2_target) * R2_target^2;
%                 b = a * (1 / R2_target - 1);
%                 R2_target_final = betarnd(a, b);
                %R2_target_final = R2_target + opt.vox_variability*randn(1);
                R2_target_final = R2_target - abs(opt.vox_variability*randn(1));
            end
        end
        
        [Y, B, R2, ~] = createSingVoxBehav( X, R2_target_final, ...
            opt.R2_tollerance, opt.vx_std_range, logger);
        voxels( :, i) = Y;
        %betas( :, i) = B;
        R2s(i) = R2;
        %pv_param(i) = pv;
    end
    if opt.vox_variability > 0 % if there is variance within subject (fpr clusterbased)
        sorted_indices = cubic_radial_sort(nb_voxels);
        [R2s, outdata_indices] = sort(R2s, "descend");
        R2s = R2s(sorted_indices); %reshape(R2s, cube_size(1), cube_size(1), cube_size(1))
        %pv_param = pv_param(outdata_indices(sorted_indices)); % !!! DOUBLE CHECK INDICES OF INDICES !!!
        %betas = betas(:, outdata_indices(sorted_indices));
        voxels = voxels(:, outdata_indices(sorted_indices));
    end
end

