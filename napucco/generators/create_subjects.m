function [p_values, R2s] = create_subjects( ...
    nb_subj, hyparams, dist_settings, ...
    saving, logger, greater_dist, save_res)
    addpath('./utils')
    addpath('./generators')
    p_values=nan; R2s= nan;

    folder_eff = hyparams.folders.SUBJDIREFF;
    folder_null = hyparams.folders.SUBJDIRNULL;
    if folder_eff ~= "none" && ~isfolder(folder_eff); mkdir(folder_eff); end
    if folder_null ~= "none" && ~isfolder(folder_null); mkdir(folder_null); end

    % create subjects data and H0es and obtain a p-value matrix from them
    if save_res
        R2s.eff = nan(nb_subj, hyparams.nb_outfeat);
        R2h0es.eff = nan(nb_subj, 3);
        p_values.eff = nan(nb_subj, hyparams.numb_permutations + 1, ...
            hyparams.nb_outfeat);
        R2s.null = nan(nb_subj, hyparams.second_cluster);
        R2h0es.null = nan(nb_subj, 3);
        p_values.null = nan(nb_subj, hyparams.numb_permutations + 1, ...
             hyparams.second_cluster);
    end
    tic;
    parfor subj_n = 1 : nb_subj
        %if dir is empty this still works as it returns empty data
        pattern = sprintf('subj_0*%d[^0-9]', subj_n);
        filelist = dir(folder_eff); %if set to none it returns empty struct!
        filelist = filelist(cellfun(@(x) ~isempty(regexp(x, pattern)), {filelist.name}));
        loading_file = fullfile(folder_eff, filelist.name);

        filelist_null = dir(folder_null); %if set to none it returns empty struct!
        filelist_null = filelist_null(cellfun(@(x) ~isempty(regexp(x, pattern)), {filelist_null.name}));
        loading_file_null = fullfile(folder_null, filelist_null.name);

        hathat = loading_file;
        if loading_file == "none"
            hathat = loading_file_null;
        end

        if isfile(hathat) % if filelist failes loading_file it is a dir
           continue
        end
        %logger.clear()
        %logger.printcr_wait("... creating subject n %d/%d", subj_n, nb_subj);

        %generate subject
        R2_target = getR2_fromdist( dist_settings.rsquared_target, greater_dist);
        [peff, pnull, ~, data , datanull] = subject_creator( ...
            hyparams, dist_settings, R2_target, folder_eff, folder_null, ...
            compose("subj_%04d", subj_n), greater_dist, subj_n, saving, logger);
%         if save_res
%             p_values.null(subj_n, :, :) = pnull;
%             p_values.eff(subj_n, :, :) = peff;
%             R2s.eff(subj_n , :) = data.R2s;
%             R2s.null(subj_n , :) = datanull.R2s;
%         end
    end
    logger.clear()
    logger.println("created %d subjects", nb_subj);
end

function R2_target = getR2_fromdist(rsquared_target, greater_dist )
    R2_target = 0;
    while ~isempty(rsquared_target) && R2_target<=0.005 || R2_target>0.6
        if greater_dist && unifrnd(0,1)>0.3
            R2_target = pearsrnd( ...
                rsquared_target(1), ...
                rsquared_target(2)*3, ...
                0.4, 2.7); %k,s: k>s**2+1
        elseif greater_dist
            R2_target = unifrnd(rsquared_target(1), ...
                rsquared_target(1) + ...
                rsquared_target(2)*5);
        else
            R2_target = rsquared_target(1) + rsquared_target(2)*randn(1);
        end
    end
end

function [param_p_values, param_p_valuesNull, subj_stats, data, datanull] = subject_creator( ...
    hyparams, dist_settings, rsquared_target, folderID, folderIDnull, ...
    subj_id, greater_dist, subj_n, saving, logger)

    subj_stats = nan;

    %generate behavioural data
    [X, ~] = create_behav(hyparams.nb_infeat, hyparams.nb_timesteps,...
            dist_settings.subj_mean_range, dist_settings.subj_std_range);

    %from behavioural data generate effected voxels activity (first cluster)
    [voxels, betas, R2s, ~] = create_cluster(X, ...
            rsquared_target, ...
            hyparams.nb_outfeat, logger, ...
            'R2_tollerance', dist_settings.threshold, ...
            'rnd_is_normal', 1, ...
            'vx_std_range', dist_settings.vox_std_range, ...
            'vox_variability', dist_settings.vox_variability);

    %if specified, generate a second cluster of voxels with no effect
    [voxels_null, R2s_null, ~] = create_nullcluster( ...
        X, dist_settings, hyparams);

    %get the null hypothesis by shuffling timepoints
    [~, R2s_H0] = get_h0es(voxels, X, hyparams.shuffler, logger);
    [~, R2s_nullH0] = get_h0es(voxels_null, X, hyparams.shufflernull, logger);
    %get p_values for the real and suffled data
%     p_values = get_pvalues(R2s, R2s_H0);
%     p_valuesNull = get_pvalues(R2s_null, R2s_nullH0);
    param_p_values = get_param_pvalues(R2s, R2s_H0, hyparams.nb_infeat, hyparams.nb_timesteps);
    param_p_valuesNull = get_param_pvalues(R2s_null, R2s_nullH0, hyparams.nb_infeat, hyparams.nb_timesteps);

    data = struct( 'subject', struct("voxels", voxels, "X",  X, "betas", betas), ...
        'p_values', param_p_values, 'R2s', R2s);
    datanull = struct( 'subject', struct("voxels", voxels_null, "X",  X, "betas", betas), ...
        'p_values', param_p_valuesNull, 'R2s', R2s_null);

    if saving && ~isempty(voxels)
        save( compose("%s.mat", fullfile(folderID, subj_id)), 'data');
        fileID = fopen( fullfile(folderID, "stats.txt"), 'a');
        fprintf(fileID, '%s %.6f\n', subj_id, mean(R2s) );
        fclose(fileID);
        fileID = fopen( fullfile(folderID, "histR2.csv"), 'a');
        fprintf(fileID, '%d, %.8f, %.8f, %.8f, %.8f\n', subj_n, mean(R2s), std(R2s), min(R2s) , max(R2s) );
        fclose(fileID);
    end
    if saving && ~isempty(voxels_null)
        save( compose("%s.mat", fullfile(folderIDnull, subj_id)), 'datanull');
        fileID = fopen( fullfile(folderIDnull, "histR2null.csv"), 'a');
        fprintf(fileID, '%d, %.8f, %.8f, %.8f, %.8f\n', subj_n, mean(R2s_null), std(R2s_null), min(R2s_null), max(R2s_null));
        fclose(fileID);
        %subj_stats = cat(1, R2s, R2s_H0);
    end
end

function [B, R2s_h0es] = get_h0es(voxels, X, shuffler, logger)
    B = nan( size(shuffler, 1), size(voxels, 2) );
    R2s_h0es = nan( size(shuffler, 1), size(voxels, 2) );
    %betas_h0es = nan( size(voxels, 2), size(shuffler, 1), size(X, 2) );
    %nvoxels = size(voxels, 2);
    if size(voxels, 2) == 0
        return
    end
    for j = 1 : size(shuffler, 1) %across the numb of perm
        %logger.print("; getting h0: %d/%d\r", j, size(shuffler, 1));
        y = voxels( shuffler(j, :), :);

        B = X\y;
        yhat = X*B;
        ssres = sum((y-yhat).^2);
        sstot = sum((y-mean(y)).^2);

        R2s_h0es(j,:) = 1-(ssres./sstot);

% 	    for vox_n = 1: nvoxels
%             [betas_h0es(vox_n, j, :), ~, ~, ~, STATS] = regress(voxels( shuffler(j, :), vox_n), X);
%             R2s_h0es(j, vox_n) = STATS(1);
%         end
    end
end

function p_values = get_pvalues(R2s, R2s_h0es)
    allR2 = cat(1, R2s, R2s_h0es);
    p_values = tiedrank(-allR2) / (size(allR2, 1) ); %maybe omit R2s row during ranking
end

function uvalue = get_param_pvalues(R2s, R2s_h0es, n_pred, n_tps)
    allR2 = cat(1, R2s, R2s_h0es);
    uvalue = 1-betacdf(allR2, double(n_pred/2), double(n_tps-n_pred-1)/2);
end