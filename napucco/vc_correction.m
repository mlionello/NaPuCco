function [hyparams, results, datadist, eff_loader, null_loader] =  vc_correction(correction, nb_subj, nb_infeat, nb_outfeat, nb_timesteps, rsquared_target, opts)
% Function Name: vc_correction
%
% Description:
%   This function is used to create subjects and to simulate experiments 
%   while performing correction methods (voxelwise, cluster-based, or
%   both). Simulate experiments should be run by interating this function N
%   times (see subj_subsampling.m)
%
% Syntax:
%   [hyparams, results, datadist, eff_loader, null_loader] = vc_correction(correction, nb_subj, nb_infeat, nb_outfeat, nb_timesteps, rsquared_target, opts)
%
% Input Arguments:
%   1. correction - String specifying the correction type ('voxelwise', 'clusterbased', 'none', 'both').
%   2. nb_subj - Number of subjects.
%   3. nb_infeat - Number of regressors.
%   4. nb_outfeat - Number of output features (must satisfy cubic root condition for cluster correction).
%   5. nb_timesteps - Number of time steps.
%   6. rsquared_target - Target R-squared value for data generation.
%   7. opts - Structure containing various optional hyperparameters and options (see the function for details).
%
% Output Arguments:
%   1. hyparams - Structure containing hyperparameters and folder information.
%   2. results - Results of the correction method.
%   3. datadist - Structure containing data distribution information.
%   4. eff_loader - Structure for effective loader.
%   5. null_loader - Structure for null loader.
%
%   Example:
%   Creation of 120 volumes with 27 effected and 80K non-effected voxels.
%   The volumes are separately saved in
%   data/in006_t0120/{nulldstXX/,r2_070_020_XX/}.
%
%   correction_type = 'none';
%   num_subjects = 120;
%   num_input_features = 6;
%   num_output_features = 27;
%   num_time_steps = 120;
%   target_rsquared = [0.07, 0.02];
%   options = struct('sec_cluster_voxels', 80000, 'numb_permutations', 2000, 'saving', 1);
% 
%   [hyparams, results, datadist, eff_loader, null_loader] = vc_correction(correction_type, num_subjects, num_input_features, num_output_features, num_time_steps, target_rsquared, options);
%
%   When simulating experiments please use subj_subsampling.m

arguments
        % correction type
        correction string {mustBeMember(correction, {'voxelwise', 'clusterbased', 'none', 'both'})} = 'none';
        
        % hyperparameters
        nb_subj int32 = 20;
        nb_infeat int32 = nan; 
        nb_outfeat int32 {musthaveintcubicroot(correction, nb_outfeat)} = 0; 
        nb_timesteps int32 = nan;
        
        % distribution values
        rsquared_target double = nan;

        % optionals hyperparams
        opts.numb_permutations int32 = 2000;
        opts.sec_cluster_voxels int32 = 0;
        
        % optional distr
        opts.threshold double = 0.002;
        opts.subj_mean_range double = [-2, 2];
        opts.subj_std_range double = [1, 2];
        opts.vox_std_range double = 1;
        opts.vox_variability = 0;
        opts.pies_thr = [0.05 0.01 0.001];
        opts.greater_dist = 0;

        opts.subj_toload = [];
        opts.saving = 1;
	    opts.save_res = 0;
        opts.eff_range_vx = [];
        opts.noneff_range_vx = [];
        opts.noneff_range_clus = [];
        opts.padding_range = [];
        opts.subj_range = [];
        opts.fpr = 1;
        opts.power = 1;

        opts.prev_settings = [];
        opts.logger = nan;
        opts.results_fold = nan;
        opts.exp_id = nan;

        % Clusterbased correction data
        opts.eff_loader = struct('init', false);
        opts.null_loader = struct('init', false);
        opts.K_eff_range = [];
        opts.K_padding_range = [];
        opts.stdev_sample_fromnull = nan;
        opts.cdt = 0.001;
    end
    addpath('./utils/')
    addpath('./generators/')
    addpath('./fwc/')

    % Initializing the output variables
    results = nan; datadist = struct();
    eff_loader = opts.eff_loader; null_loader = opts.null_loader;

    % Setting up the variables and the folders
    pies_thr = opts.pies_thr;
    current_run_id = string( datetime( 'now', 'Format', 'yyMMdd_HHmmss'));
    if isempty(opts.prev_settings)
        [hyparams, dist_settings] = prepeare_data( ...
                opts, nb_outfeat, nb_subj, nb_timesteps, ...
                nb_infeat, rsquared_target, opts.greater_dist);
    else
        hyparams = opts.prev_settings.hyparams;
        dist_settings = opts.prev_settings.dist_settings;
    end

    % Creating a logger
    logfile = fullfile(hyparams.folders.LOGS, "log_" + current_run_id + ".txt");
    if isnan(opts.logger)
        logger = MPrint(logfile);
    else
        logger = opts.logger;
    end
    logger.println("working directory: %s", hyparams.folders.HOMEPATH);

    % Generating nb_subj subjects
    if ~isempty(opts.subj_toload)
        [p_values, R2s, eff_loader, null_loader] = subj_loader(hyparams, logger, opts.subj_toload, eff_loader, opts.null_loader);
    else
        [p_values, ~] = create_subjects(nb_subj, ...
            hyparams, dist_settings, opts.saving, logger, ...
            opts.greater_dist, opts.save_res);
    end

    % If no correction are applied (only subject data generation), exit!
    if correction == "none"; return; end

    if ~isnan(opts.stdev_sample_fromnull) % HOW TO OPTIMIZE?????! 
        resamping_ind = resample_nullR2(R2s.null, opts.stdev_sample_fromnull, opts.noneff_range_vx(end));
        parfor subj_i = 1: nb_subj
            resampledr2(subj_i, :) = R2s.null(subj_i, resamping_ind(subj_i, :));
            resampledp(subj_i, :, :) = p_values.null(subj_i, :, resamping_ind(subj_i, :));
        end
        R2s.null = resampledr2;
        p_values.null = resampledp;
    end

    datadist.r2 = R2s;
    datadist.p_values = p_values;

    % Non-parametric combination
    logger.print("COMBINING p_values... ");
    if isempty(opts.subj_range)
        opts.subj_range = 5:5:hyparams.nb_subj;
    end
    superman = nonpc(p_values, hyparams, opts.subj_range );
    datadist.superman = superman;
    datadist.subj_range = opts.subj_range;
    logger.println("Done -> superman");

    logger.println("-- STARTING CORRECTION METHOD --");
    if ~isfolder(hyparams.folders.RESULTS); mkdir(hyparams.folders.RESULTS); end

    if correction == "voxelwise" || correction == "both"
        % Voxelwise correction
        logger.println("VOXELWISE correction: ");
        results = voxelwise_correction(superman, ...
            hyparams, opts.eff_range_vx, opts.noneff_range_vx, opts.subj_range, pies_thr, opts, logger);
        if opts.save_res
            outfolder = fullfile( opts.results_fold, "voxelwise_correction", 'grouping');
            if ~exist(outfolder, 'dir'); mkdir(outfolder); end
		    save(fullfile( outfolder, opts.exp_id), 'results');
		    %save(fullfile( results_fold, "datadist"), 'datadist');
        end 
    end

    if correction == "clusterbased" || correction == "both"
        % Clusterbased correction
        logger.println("CLUSTER BASED correction... ");
        results = cluster_correction(superman, ...
            opts.subj_range, pies_thr, opts, logger, opts.cdt, ...
            opts.noneff_range_clus, opts.padding_range);
        if opts.save_res
            outfolder = fullfile( opts.results_fold, "clusterbased_correction", 'grouping');
            if ~exist(outfolder, 'dir'); mkdir(outfolder); end
		    save(fullfile( outfolder, opts.exp_id), 'results');
		    %save(fullfile( results_fold, "datadist"), 'datadist');
        end
    end

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUPPORT FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hyparams, dist_settings] = prepeare_data( ...
                opts, nb_outfeat, nb_subj, nb_timesteps, nb_infeat, rsquared_target, grt_dist)
    %setting fodler structure
    folders = conf_folders( nb_infeat, nb_timesteps, ...
        rsquared_target, grt_dist, opts.vox_variability, [nb_outfeat, opts.sec_cluster_voxels]);

    [hyparams, dist_settings] = setparameters( ...
        nb_outfeat, opts, nb_subj, nb_timesteps, ...
        nb_infeat, rsquared_target, folders);

    % saving settings in the subject folders to be used
    subjdir = [ hyparams.folders.SUBJDIREFF, hyparams.folders.SUBJDIRNULL];
    for fold_i = 1:2
        if subjdir(fold_i) ~= "none"
        save( fullfile( './', subjdir(fold_i), 'hyparams'), 'hyparams');
        save( fullfile( './', subjdir(fold_i), 'dist_settings'), 'dist_settings');
        end
    end
end

function musthaveintcubicroot(correction, a)
    % Test for cubic shape of clusters - needed for cluster correction.
    if correction ~= "clusterbased"
        return
    end
    tmp = nan(a, 1);
    cubirootdim = int32(power(double(a), 1/3));
    try
        tmp = reshape(tmp, cubirootdim, cubirootdim, cubirootdim);
    catch
        eid = 'Size:notACube';
        msg = 'The cubic root of the size must be an integer';
        throwAsCaller(MException(eid,msg))
    end
end

