function hyparams = generate_volumes(nb_volumes, nb_infeat, nb_timesteps, nb_outfeat, rsquared_target, opts)
% Function Name: generate_volumes
%
% Description:
%   This function is used to create volumes ...
%
% Syntax:
%   generate_volumes( nb_volumes, nb_infeat, nb_timesteps, nb_outfeat, rsquared_target, opts)
%
% Input Arguments:
%   1. nb_volumes - Number of volumes.
%   2. nb_infeat - Number of regressors.
%   3. nb_timesteps - Number of time steps.
%   4. nb_outfeat - Number of output features (must satisfy cubic root condition for cluster correction).
%   5. rsquared_target - Target R-squared value for data generation.
%   6. opts - Structure containing various optional hyperparameters and options (see the function for details).
%
%   Example:
%   Creation of 120 volumes with 27 effected and 80K non-effected voxels.
%   The volumes are separately saved in
%   data/in006_t0120/{noneffdstXX/,r2_070_020_XX/}.
%
%   nb_volumes = 120;
%   num_input_features = 6;
%   num_output_features = 27;
%   num_time_steps = 120;
%   target_rsquared = [0.07, 0.02];
%   options = struct('nb_noneffvx', 80000, 'numb_permutations', 2000, 'saving', 1);
% 
%   generate_volumes(nb_volumes, num_input_features, num_time_steps, num_output_features, target_rsquared, options);

arguments
        nb_volumes int32;
        nb_infeat int32;
        nb_timesteps int32;
        
        nb_outfeat int32 = 0;
        rsquared_target double = nan;

        % optionals hyperparams
        opts.numb_permutations int32 = 1000;
        opts.nb_noneffvx int32 = 0;
        
        % optional distr
        opts.threshold double = 0.002;
        opts.subj_mean_range double = [-2, 2];
        opts.subj_std_range double = [1, 2];
        opts.vox_std_range double = 1;
        opts.vox_variability = 0;
        opts.greater_dist = 0;

        opts.prev_settings = [];
        opts.logger = nan;
        opts.exp_id = nan;
        opts.rootdir = './';
    end
    addpath('./utils/')
    addpath('./generators/')
    addpath('./fwc/')

    % Setting up the variables and the folders
    current_run_id = string( datetime( 'now', 'Format', 'yyMMdd_HHmmss'));
    if isempty(opts.prev_settings)
        [hyparams, dist_settings] = prepeare_data( ...
                opts, nb_outfeat, nb_volumes, nb_timesteps, ...
                nb_infeat, rsquared_target, opts.greater_dist, opts.rootdir);
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

    folder_eff = hyparams.folders.SUBJDIREFF;
    folder_noneff = hyparams.folders.SUBJDIRNONEFF;
    if folder_eff ~= "none" && ~isfolder(folder_eff); mkdir(folder_eff); end
    if folder_noneff ~= "none" && ~isfolder(folder_noneff); mkdir(folder_noneff); end

    tic;
    parfor_progress(nb_volumes);
    parfor vol_n = 1 : nb_volumes
        parfor_progress;
        if check_file_already_in(vol_n, folder_eff, folder_noneff)
            continue
        end

        % sample R2 target
        R2_target = getR2_fromdist( dist_settings.rsquared_target, opts.greater_dist);

        % generate volume
        subject_creator( ...
            hyparams, dist_settings, R2_target, folder_eff, folder_noneff, ...
            compose("vol_%04d", vol_n), vol_n, logger);
    end
    parfor_progress(0);
    logger.clear()
    logger.println("created %d volumes", nb_volumes);
end

function R2_target = getR2_fromdist(rsquared_target, greater_dist )
    R2_target = 0;
    while ~isempty(rsquared_target) && R2_target<=0.005 || R2_target>0.6
        if greater_dist && unifrnd(0,1)>0.3
            R2_target = pearsrnd( ...
                rsquared_target(1), ...
                rsquared_target(2)*3, ...
                0.4, 2.7);
        elseif greater_dist
            R2_target = unifrnd(rsquared_target(1), ...
                rsquared_target(1) + ...
                rsquared_target(2)*5);
        else
            R2_target = rsquared_target(1) + rsquared_target(2)*randn(1);
        end
    end
end

function [param_p_values, param_p_noneffvalues, subj_stats, data_eff, data_noneff] = subject_creator( ...
    hyparams, dist_settings, rsquared_target, folderID, folderIDnoneff, ...
    vol_id, vol_n, logger)

    subj_stats = nan;

    % generate behavioural data
    [X, ~] = create_behav(hyparams.nb_infeat, hyparams.nb_timesteps,...
            dist_settings.subj_mean_range, dist_settings.subj_std_range);

    % from behavioural data generate effected voxels activity (first cluster)
    [voxels_eff, betas, R2s, ~] = create_cluster(X, ...
            rsquared_target, ...
            hyparams.nb_outfeat, logger, ...
            'R2_tollerance', dist_settings.threshold, ...
            'rnd_is_normal', 1, ...
            'vx_std_range', dist_settings.vox_std_range, ...
            'vox_variability', dist_settings.vox_variability);

    % if specified, generate a second cluster of voxels with no effect
    [voxels_noneff, R2s_noneff, ~] = create_noneffcluster( ...
        X, dist_settings, hyparams);

    % get the null hypothesis by shuffling timepoints
    [~, R2s_H0] = get_h0es(voxels_eff, X, hyparams.shuffler);
    [~, R2s_noneffH0] = get_h0es(voxels_noneff, X, hyparams.shufflernoneff);
    % get p_values for the real and suffled data
    param_p_values = get_param_pvalues(R2s, R2s_H0, hyparams.nb_infeat, hyparams.nb_timesteps);
    param_p_noneffvalues = get_param_pvalues(R2s_noneff, R2s_noneffH0, hyparams.nb_infeat, hyparams.nb_timesteps);

    data_eff = struct( 'subject', struct("voxels", voxels_eff, "X",  X, "betas", betas), ...
        'plog_values', log(single(param_p_values)));
    data_noneff = struct( 'subject', struct("voxels", voxels_noneff, "X",  X, "betas", betas), ...
        'plog_values', log(single(param_p_noneffvalues)));

    if ~isempty(voxels_eff)
        save( compose("%s.mat", fullfile(folderID, vol_id)), 'data_eff');
        fileID = fopen( fullfile(folderID, "stats.txt"), 'a');
        fprintf(fileID, '%s %.6f\n', vol_id, mean(R2s) );
        fclose(fileID);
        fileID = fopen( fullfile(folderID, "histR2.csv"), 'a');
        fprintf(fileID, '%d, %.8f, %.8f, %.8f, %.8f\n', vol_n, mean(R2s), std(R2s), min(R2s) , max(R2s) );
        fclose(fileID);
    end
    if ~isempty(voxels_noneff)
        save( compose("%s.mat", fullfile(folderIDnoneff, vol_id)), 'data_noneff');
        fileID = fopen( fullfile(folderIDnoneff, "histR2noneff.csv"), 'a');
        fprintf(fileID, '%d, %.8f, %.8f, %.8f, %.8f\n', vol_n, mean(R2s_noneff), std(R2s_noneff), min(R2s_noneff), max(R2s_noneff));
        fclose(fileID);
    end
end

function [B, R2s_h0es] = get_h0es(voxels, X, shuffler)
    B = nan( size(shuffler, 1), size(voxels, 2) );
    R2s_h0es = nan( size(shuffler, 1), size(voxels, 2) );
    if size(voxels, 2) == 0
        return
    end
    for j = 1 : size(shuffler, 1) % across the numb of perm
        y = voxels(shuffler(j, :), :);

        B = X \ y;
        yhat = X * B;
        ssres = sum((y - yhat) .^ 2);
        sstot = sum((y - mean(y)) .^ 2);

        R2s_h0es(j, :) = 1 - (ssres ./ sstot);

    end
end

function uvalue = get_param_pvalues(R2s, R2s_h0es, n_pred, n_tps)
    allR2 = cat(1, R2s, R2s_h0es);
    uvalue = 1-betacdf(allR2, double(n_pred/2), double(n_tps-n_pred-1)/2);
end

function is_in = check_file_already_in(vol_n, folder_eff, folder_noneff)
        pattern = sprintf('vol_0*%d[^0-9]', vol_n);
        filelist = dir(folder_eff); % if set to none it returns empty struct!
        filelist = filelist(cellfun(@(x) ~isempty(regexp(x, pattern)), {filelist.name}));
        loading_file = fullfile(folder_eff, filelist.name);

        filelist_noneff = dir(folder_noneff); % if set to none it returns empty struct!
        filelist_noneff = filelist_noneff(cellfun(@(x) ~isempty(regexp(x, pattern)), {filelist_noneff.name}));
        loading_file_noneff = fullfile(folder_noneff, filelist_noneff.name);

        is_in = loading_file == "none" && isfile(loading_file_noneff);
end


function [hyparams, dist_settings] = prepeare_data( ...
                opts, nb_outfeat, nb_subj, nb_timesteps, nb_infeat, rsquared_target, grt_dist, rootdir)
    % setting fodler structure
    folders = conf_folders( nb_infeat, nb_timesteps, ...
        rsquared_target, grt_dist, opts.vox_variability, [nb_outfeat, opts.nb_noneffvx], rootdir);

    [hyparams, dist_settings] = setparameters( ...
        nb_outfeat, opts, nb_subj, nb_timesteps, ...
        nb_infeat, rsquared_target, folders);

    % saving settings in the subject folders to be used
    subjdir = [ hyparams.folders.SUBJDIREFF, hyparams.folders.SUBJDIRNONEFF];
    for fold_i = 1 : 2
        if subjdir(fold_i) ~= "none"
        save( fullfile( subjdir(fold_i), 'hyparams'), 'hyparams');
        save( fullfile( subjdir(fold_i), 'dist_settings'), 'dist_settings');
        end
    end
end