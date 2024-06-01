function [hyparams, results, ...
    eff_loader, noneff_loader] =  vc_correction(correction, ...
    nb_subj, subj_toload, eff_loader, noneff_loader, settings, opts)
% Function Name: vc_correction
%
% Description:
%   This function to simulate a single experiment
%   by performing correction methods (voxelwise, cluster-based, or
%   both). Simulate experiments should be run by interating this function N
%   times (see subj_subsampling.m)
%
% Syntax:
%   [hyparams, results, eff_loader, noneff_loader] = vc_correction(correction, nb_subj, nb_infeat, nb_outfeat, nb_timesteps, rsquared_target, opts)
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
%   4. eff_loader - Structure for effective loader.
%   5. noneff_loader - Structure for noneff loader.
%
%   Calling directly this method is deprecated; when simulating experiments, please use subj_subsampling.m

arguments
    % correction type
    correction string {mustBeMember(correction, {'voxelwise', 'clusterbased', 'both'})} = 'both';
    nb_subj = 100;
    subj_toload = [];

    eff_loader = struct('init', false);
    noneff_loader = struct('init', false);

    settings = [];

    opts.save_res = 0;
    opts.eff_range_vx = [];
    opts.noneff_range_vx = [];
    opts.noneff_range_clus = [];
    opts.padding_range = [];
    opts.subj_range = [];
    opts.fpr = 1;
    opts.power = 1;
    opts.pies_thr = [0.05, 0.01, 0.001];

    opts.logger = nan;
    opts.results_fold = nan;
    opts.exp_id = nan;

    % Clusterbased correction data
    opts.K_eff_range = [];
    opts.K_padding_range = [];
    opts.stdev_sample_fromnoneff = nan;
    opts.cdt = 0.001;
    opts.p_noneff = [0];
    opts.n_additional_row = 0;
    opts.rootdir = './';
    opts.cpus = 1;
end
    addpath('./utils/')
    addpath('./generators/')
    addpath('./fwc/')

    % Initializing the output variables
    results = nan;

    % Setting up the variables and the folders
    pies_thr = opts.pies_thr;
    current_run_id = string( datetime( 'now', 'Format', 'yyMMdd_HHmmss'));

    hyparams = settings.hyparams;
    dist_settings = settings.dist_settings;

    % Creating a logger
    logfile = fullfile(hyparams.folders.LOGS, "log_" + current_run_id + ".txt");
    if isnan(opts.logger)
        logger = MPrint(logfile);
    else
        logger = opts.logger;
    end
    logger.println("working directory: %s", hyparams.folders.HOMEPATH);

    [p_values, ~, eff_loader, noneff_loader] = subj_loader(hyparams, logger, subj_toload, eff_loader, noneff_loader);

    % Non-parametric combination
    logger.print("COMBINING p_values... ");
    if isempty(opts.subj_range)
        opts.subj_range = 5:5:hyparams.nb_subj;
    end
    superman = nonpc(p_values, hyparams, opts.subj_range, opts.n_additional_row, opts.p_noneff, opts.cpus);
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
        end
    end

    if correction == "clusterbased" || correction == "both"
        % Clusterbased correction
        logger.println("CLUSTER BASED correction... ");
        results = cluster_correction(superman, ...
            opts.subj_range, pies_thr, opts, logger, opts.cdt, length(opts.p_noneff), ...
            opts.noneff_range_clus, opts.padding_range);
        if opts.save_res
            outfolder = fullfile( opts.results_fold, "clusterbased_correction", 'grouping');
            if ~exist(outfolder, 'dir'); mkdir(outfolder); end
		    save(fullfile( outfolder, opts.exp_id), 'results');
        end
    end

end
