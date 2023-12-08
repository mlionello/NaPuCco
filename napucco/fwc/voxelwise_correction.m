function results = voxelwise_correction(superman, hyparams, eff_range, non_eff_range, subj_range, pies_thr, opts, logger)
arguments
    superman;
    hyparams;
    eff_range;
    non_eff_range;
    subj_range;
    pies_thr;
    opts;
    logger;
end
    if isempty(eff_range)
        eff_range = 1;
    end

    if isempty(non_eff_range)
        non_eff_range = 500: 500: hyparams.second_cluster;
        non_eff_range = unique([non_eff_range, hyparams.second_cluster]);
    end

    results = struct();
    logger.print("obtaining FWC... ");
    fwc_mixed = nan( eff_range(end), size(superman.eff,1), length(eff_range),  length(non_eff_range));
    fwc_null = nan( non_eff_range(end), size(superman.null, 1), length(non_eff_range));

    % FPR
    if opts.fpr
        logger.print("\tgetting FPR from non-correlated data... ");
        for v = 1: length(non_eff_range)
            %logger.print('set of voxels nr %d/%d\r', v, length(non_eff_range))
            [tmp, ~, ~] = get_fwc( ...
	            superman.null(:, :, 1: non_eff_range(v))); % this returns a [nb_voxels, nb_subj] matrix
            fwc_null(1: non_eff_range(v), :, v) = tmp;
            %fwc_null( :, v) = tmp(1,:); 
        end
        results.fwc_null = fwc_null;
        for p_i = 1: length(pies_thr) % be careful with nan values
            results.fpr(1:size(fwc_null,1), ...
                1:size(fwc_null,2), ...
                1:size(fwc_null,3), p_i) = fwc_null < pies_thr(p_i);
        end
        logger.println("Done -> results.fwc_null, results.fpr");
    end

    % POWER
    if opts.power
        logger.print("\tgetting POWER from mixed data... ");
        for v_null = 1: length(non_eff_range)
            logger.print('set of voxels nr %d/%d\r', v_null, length(non_eff_range))
            for v_eff = 1: length(eff_range)
                ultimate_brain = cat(3, ...
                    superman.eff(:, :, eff_range(v_eff)), ...
                    superman.null(:, :, 1: non_eff_range(v_null)));

                [tmp , ~, ~] = get_fwc(ultimate_brain);

                % From tmp I extract the info fwc from only the
                % effected voxels, discarding the non-effected ones;
                % In the case only 1 voxel is effected, it returns a
                % matrix nb_subj times length(range of non-eff voxels).
                fwc_mixed( ...
                    1: eff_range(v_eff), :, v_eff, v_null) = tmp(1: eff_range(v_eff), :);
            end
        end
        results.fwc_mixed = fwc_mixed;
        % fwc_mixed is already tiedranked(-x). The resutls lower then
        % the threshold are passed experiments, higher the
        % threshold are failed.
        % The result from each effected voxels is compared with the
        % threhsold.
        for p_i = 1: length(pies_thr)
            results.power(1:size(fwc_mixed,1), ...
                1:size(fwc_mixed,2), ...
                1:size(fwc_mixed,3), ...
                1:size(fwc_mixed,4), p_i) = fwc_mixed < pies_thr(p_i);
        end
        logger.println("Done -> results.fwc_mixed, results.power");
    end

    results.pies_thr = pies_thr;
    results.subj_range = subj_range;
    results.eff_range = eff_range;
    results.non_eff_range = non_eff_range;
end