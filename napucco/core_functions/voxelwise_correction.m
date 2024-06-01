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
        non_eff_range = 500: 500: hyparams.nb_noneffvx;
        non_eff_range = unique([non_eff_range, hyparams.nb_noneffvx]);
    end
    maxnoneffvx = non_eff_range(end);

    results = struct();
    fwc_mixed = cell( length(non_eff_range));
    fwc_noneff = cell( length(non_eff_range), 1);
    results.fwc_noneff = nan( non_eff_range(end), size(superman.noneff, 1), length(non_eff_range));
    results.fwc_mixed = nan( size(superman.noneff, 1), length(non_eff_range) );
    results.fpr = nan( size(superman.noneff, 1), length(non_eff_range), length(pies_thr));
    results.power = nan(size(superman.noneff, 1), length(non_eff_range), length(pies_thr));
    % FPR
    noneff_superman = superman.noneff(:, :, 1:maxnoneffvx);
    if opts.fpr
        logger.print("obtaining FWC... ");
        logger.println("\tgetting FPR from non-correlated data... ");
        parfor_progress(length(non_eff_range)); tic;
        for v = 1: length(non_eff_range)
            parfor_progress;

            tmp = get_fwc( ...
	            noneff_superman(:, :, 1 : non_eff_range(v)), pies_thr); % this returns a [nb_voxels, nb_subj] matrix
            results.fpr( :, v, :) = any(tmp);

        end
        parfor_progress(0);
        toc;
        logger.println("Done -> results.fwc_noneff, results.fpr");
    end

    % POWER
    if opts.power
        logger.println("\tgetting POWER from mixed data... ");
        eff_slice = superman.eff(:, :, 1);
        [N, M, ~] = size(noneff_superman);
        parfor_progress(length(non_eff_range)); tic;
        for v_noneff = 1: length(non_eff_range)
                parfor_progress;

                current_range = non_eff_range(v_noneff);
                ultimate_brain = zeros(N, M, current_range + 1);
                ultimate_brain(:, :, 1) = eff_slice;
                ultimate_brain(:, :, 2:end) = noneff_superman(:, :, 1:current_range);

                tmp = get_fwc(ultimate_brain, pies_thr);
                results.power(:, v_noneff, :) = squeeze(tmp(1, :, :)); % keep pvalue only for the effected VX

        end
        parfor_progress(0);
        toc;
        logger.println("Done -> results.fwc_mixed, results.power");
    end

    results.pies_thr = pies_thr;
    results.subj_range = subj_range;
    results.eff_range = eff_range;
    results.non_eff_range = non_eff_range;
end