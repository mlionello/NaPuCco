function [voxels_eff_pv, voxels_null_pv] = get_chi2_dist_from_cluster(voxels_eff, voxels_null, parametric, subj_range, len_pnull)
    repeat_nu_on_eff = @(size_last_dimension) repmat(double(repmat(2 * subj_range, 1, len_pnull)), ...
                    size(voxels_eff, 1), ...
                    1, ...
                    size_last_dimension);
    repeat_nu_on_null = @(size_last_dimension) repmat(double(repmat(2 * subj_range, 1, len_pnull)), ...
                    size(voxels_null, 1), ...
                    1, ...
                    size_last_dimension);

    if parametric
        voxels_eff_pv = 1 - chi2cdf(voxels_eff, repeat_nu_on_eff(size(voxels_eff, 3)));

        max_cpus = 4;
        nb_batches = 12;
        batchsize = int32(floor(size(voxels_null, 3) / nb_batches));

        numel_indices = cell(nb_batches, 1);
        voxels_null_batch = cell(nb_batches, 1);
        for i = 1 : nb_batches
            onset = (i - 1) * batchsize + 1;
            offset = i * batchsize;
            indices = onset : offset;
            voxels_null_batch{i} = voxels_null(:, :, indices);
            numel_indices{i} = numel(indices);
        end

        tic;
        parfor_progress(nb_batches);
        parfor (i = 1:nb_batches, min([nb_batches, max_cpus, gcp("nocreate").NumWorkers]))
            parfor_progress;
            voxels_null_pv{i} = 1 - chi2cdf(...
                voxels_null_batch{i}, ...
                repeat_nu_on_null(numel_indices{i}));
        end
        parfor_progress(0);
        voxels_null_pv = cat(3, voxels_null_pv{:});
        toc;

        remaining_indices = (nb_batches + 1) * batchsize : size(voxels_null, 3);
        if ~isempty(remaining_indices)
            remaining_subbatch = 1 - chi2cdf( ...
                    voxels_null(:, :, remaining_indices), ...
                    repeat_nu_on_null(numel(remaining_indices)));
    
            voxels_null_pv(:, :, remaining_indices) = squeeze(remaining_subbatch);
        end
    else
        voxels_eff_pv = tiedrank(-voxels_eff) / size(voxels_eff, 1); % MINUS IF FISHER
        voxels_null_pv = tiedrank(-voxels_null) / size(voxels_null, 1); % MINUS IF FISHER
    end
end