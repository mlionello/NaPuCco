function resfolder = subj_subsampling(correction, datapath, nb_subj, nb_rep, r2, opts)
% SUBJ_SUBSAMPLING  Applies the correction methods to subsamples of a given dataset.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(CORRECTION, PATH, NB_SUBJ, NB_REP) It applies the given
%   correction to NB_SUBJ subjects extracted randomly from the given PATH and
%   it repeats the extraction and the correction NB_REP times (1000 if not specified).
%   the results are saved in 'grouping' folder.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(CORRECTION, PATH, NB_SUBJ, NB_REP, R2) It
%   samples the NB_SUBJ sampled from a gaussian fit to the given R2.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'p_noneffsubj', p_noneffsubj) It
%   introduce a stochastic portion of non effected participants in the
%   population
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'stack_mem', stack_mem) It allows to
%   store stack_mem number of volumes (for effected and non-effected ones)
%   in a priority queue system refraining from loading the same volumes
%   multiple times over the NB_REP repetitions.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'subj_range', subj_range) It repeats the corrections for the given intervals of number
%   of participants
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'eff_range_vx', eff_range_vx) It repeats the corrections for the given intervals of number
%   of effected voxels in voxelwise correction.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'noneff_range_vx', noneff_range_vx) It repeats the corrections for the given intervals of number
%   of NON-effected voxels in voxelwise correction.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'noneff_range_clus', noneff_range_clus) It repeats the corrections for the given intervals of number
%   of noneff_range_clus^3 NON-effected voxels in FPR for cluster correction.
%
%   OUTFOLDER = SUBJ_SUBSAMPLING(..., 'padding_range', padding_range) It
%   repeats the corrections for the given intervals of padding width
%   of NON-effected voxels in Power calcluation for cluster correction.
%
arguments
    correction string {mustBeMember(correction, {'voxelwise', 'clusterbased', 'none', 'both'})} = 'none';
    datapath = [];
    nb_subj = 80;
    nb_rep = 1000;
    r2 = []; % if r2 dist is not given random participants are extracted
    opts.p_noneffsubj = 0;
    opts.power = 1;
    opts.fpr = 1;
    opts.pies = [0.05, 0.01, 0.001];
    opts.eff_range_vx = [];
    opts.noneff_range_vx = [];
    opts.noneff_range_clus = [];
    opts.padding_range = [];
    opts.subj_range = [];
    opts.vox_variability = 0;
    opts.stdev_sample_fromnoneff = nan;
    opts.resfolder = "";
    opts.cdt = 0.001;
    opts.stack_mem = 600;
    opts.save_res = 0;
    opts.rootdir = './';
    opts.cpus = 4;
end
    addpath('./utils')
    addpath('./fwc')
    addpath('./generators')

    p = gcp('nocreate');
    if isempty(p)
        parpool(opts.cpus)
    elseif p.NumWorkers ~= opts.cpus
        delete(gcp('nocreate'))
        parpool(opts.cpus)
    end
    
    % Congruency-check among the given arguments
    argscheck(opts.power, datapath);

    % Extracting information from the dataset
    [hyparams, dist_settings, subject_list, file_list] = parse_path(datapath, nb_subj);
    if ~opts.power; hyparams.folders.SUBJDIREFF = ""; end
    
    % Folders for storing the results
    [hyparams, resfolder] = get_output_folders(hyparams, r2, datapath, opts, correction);

    if opts.power && ~matches(correction, "voxelwise")
        hyparams.cdt = opts.cdt;
    end

    % Providing a structure containing all the settings to be loaded
    prev_settings = struct('hyparams', hyparams, 'dist_settings', dist_settings);
 
    % initializing loaders (stack of files
    eff_loader = struct( 'init', true, 'len', opts.stack_mem, 'stack_indices', []);
    eff_loader.data = cell(0); eff_loader.loaded_names = cell(0);  
    noneff_loader = struct('init', true, 'len', opts.stack_mem, 'stack_indices', []);
    noneff_loader.data = cell(0); noneff_loader.loaded_names = cell(0);
    
    % Starting the NB_REP of repetitions
    init_i = length(dir(resfolder)) - 1;
    if exist(fullfile(resfolder, "clusterbased_correction", 'grouping'), 'dir')
        init_i = length(dir(fullfile(resfolder, "clusterbased_correction", 'grouping'))) - 1;
    end
    if exist(fullfile(resfolder, "voxelwise_correction", 'grouping'), 'dir')
        init_i = min(init_i, length(dir(fullfile(resfolder, "voxelwise_correction", 'grouping'))) - 1);
    end

    for i = init_i: nb_rep
        if ~isempty(r2)
            fprintf("%d repetition out of %d, nb_subjects: %d, r2 mean: %.3f, std: %.3f, subj_range: [",...
                i, nb_rep, hyparams.nb_subj, r2(1), r2(2))
            fprintf("%d ", opts.subj_range)
            fprintf("]\n")
        else
            fprintf("%d repetition out of %d, nb_subjects: %d, subj_range: [",...
                i, nb_rep, hyparams.nb_subj)
            fprintf("%d ", opts.subj_range)
            fprintf("]\n")
        end

        % Loading list to be passed to create_subj, if r2 is given it
        % extracts from the primary dataset the nb_subj subjects
        % whose recorded r2 is the closest to the nb_subj sampling of
        % a gaussian fitting the r2 params.
        if ~isempty(r2)
            r2_targets = r2(1) + r2(2)*randn(hyparams.nb_subj, 1);
            [loading_list, n_additional_row, ~] = get_closest_samples( subject_list, ...
                 r2_targets, file_list, opts.p_noneffsubj, hyparams);
        else
            [loading_list, n_additional_row] = get_random_samples( opts.p_noneffsubj, hyparams);
        end
        if ~isfile(fullfile(resfolder, 'hyparams.mat'))
            save(fullfile(resfolder,'hyparams.mat'), 'hyparams')
            save(fullfile(resfolder,'opts.mat'), 'opts')
        end

        % Starting the non-parametric combination and the corrections
        initTime = tic;
        [~, results, eff_loader, noneff_loader] = vc_correction( ...
                correction, ...
                hyparams.nb_subj, ...
                loading_list, ...
                eff_loader, ...
                noneff_loader, ...
                prev_settings, ...
	            'save_res', opts.save_res, ... %saving res are done here below
                'subj_range',opts.subj_range, ...
                'eff_range_vx', opts.eff_range_vx, ...
                'noneff_range_vx', opts.noneff_range_vx, ...
                'noneff_range_clus', opts.noneff_range_clus, ...
                'padding_range', opts.padding_range, ...
                'stdev_sample_fromnoneff', opts.stdev_sample_fromnoneff, ...
                'pies_thr', opts.pies, ...
                'power', opts.power, 'fpr',opts.fpr, ...
                'results_fold', resfolder, ...
                'exp_id', compose("results_%04d.mat", i), ...
                'p_noneff', opts.p_noneffsubj, ...
                'n_additional_row', n_additional_row, ...
                'cdt', opts.cdt, ...
                'cpus', opts.cpus);
        fprintf("total time: %s\n", duration(0,0,toc(initTime)));

        % The results are saved in outfolder
        if ~opts.save_res
            if ~exist(fullfile(resfolder, 'grouping'), 'dir')
                mkdir(fullfile(resfolder, 'grouping'));
            end
            save(fullfile(resfolder, 'grouping', ...
                compose("results_%04d.mat", i)), "results");
        end

    end
end

function argscheck(power, datapath)
    if power && length(datapath) == 1
        eid = 'argscombo:notsupported';
        msg = ['if power is demanded, a ' ...
            'second path to the non effected dataset is required'];
        throwAsCaller(MException(eid,msg))
    end
    if ~power && length(datapath) > 1
        eid = 'argscombo:notsupported';
        msg = 'you set power off: a second path is given, but it was not required!';
        throwAsCaller(MException(eid,msg))
    end
end

function [hyparams, dist_settings, subject_list, file_list] = parse_path(datapath, nb_subj)
    homepath = datapath(1);

    file_list = dir(homepath);
    file_list = {file_list(cellfun(@(x) contains(x,"vol_"), {file_list.name})).name};

    subject_list = nan;    
    histr2file = fullfile(homepath, "histR2.csv");
    if isfile(histr2file)
        subject_list = csvread(histr2file);
    end

    dist_settings = load(fullfile(homepath, "dist_settings.mat"));
    dist_settings = dist_settings.dist_settings;
    load(fullfile(homepath, "hyparams.mat"), 'hyparams');
    hyparams.nb_subj = nb_subj;

    % If an additional path is given, it applies the path and the number of
    % voxels to the settings that will be passed to the vc_correction
    if length(datapath)>1
        noneffpath = datapath(2);
        hyparamsnoneff = load(fullfile(noneffpath, "hyparams.mat"));
        hyparams.folders.SUBJDIRNONEFF = hyparamsnoneff.hyparams.folders.SUBJDIRNONEFF;
        hyparams.nb_noneffvx = hyparamsnoneff.hyparams.nb_noneffvx;
        hyparams.shufflernoneff = hyparamsnoneff.hyparams.shufflernoneff;
        if hyparamsnoneff.hyparams.numb_permutations ~= hyparams.numb_permutations
            % eid = 'perm:notmatching';
            % msg = compose("the two datasets have different number of permutations: %d vs %d !", ...
            %     hyparams.numb_permutations, hyparamsnoneff.hyparams.numb_permutations);
            %             throwAsCaller(MException(eid, msg))
            hyparams.numb_permutations = min(hyparams.numb_permutations,hyparamsnoneff.hyparams.numb_permutations);
        end
    end

end


function [loading_files, n_additional_row] = get_random_samples( p_noneffsubj, hyparams)
    n_additional_row = 0;
    loading_files(:, [1,2]) = repelem("none", hyparams.nb_subj, 2);

    eff_files = dir(hyparams.folders.SUBJDIREFF);
    eff_files = string({eff_files(cellfun(@(x) contains(x, "vol_"), {eff_files.name})).name});
    if ~isempty(eff_files)
        eff_ind = randsample(length(eff_files), hyparams.nb_subj, true); % WITH replacement
        loading_files(:, 1) = fullfile(hyparams.folders.SUBJDIREFF, eff_files(eff_ind))';
    end

    noneff_files = dir(hyparams.folders.SUBJDIRNONEFF);
    noneff_files = string({noneff_files(cellfun(@(x) contains(x, "vol_"), {noneff_files.name})).name});
    if ~isempty(noneff_files)
        noneff_ind = randsample(length(noneff_files), hyparams.nb_subj, true); % WITH replacement
        loading_files(:, 2) = fullfile(hyparams.folders.SUBJDIRNONEFF, noneff_files(noneff_ind))';
    end

    if any(p_noneffsubj > 0) % noneff effected subjects
        % noneff_indices = rand(length(loading_files), 1) < max(p_noneffsubj);
        n_additional_row = int32(hyparams.nb_subj*max(p_noneffsubj));
        noneff_ind = randsample(length(noneff_files), n_additional_row);
        loading_files(hyparams.nb_subj + 1 : hyparams.nb_subj + n_additional_row, :) = [ ...
            fullfile(hyparams.folders.SUBJDIRNONEFF, noneff_files(noneff_ind))',...
            loading_files(1 : n_additional_row, 2) ...
            ];
    end
end


function [loading_files, n_additional_row, r2_real] = get_closest_samples(r2_available, r2_targets, filelist, p_noneffsubj, hyparams)
    subj_indices = [];
    n_additional_row = 0;
    for r2_i = 1: length(r2_targets)
        [~, subj_indices(r2_i)] = min(abs(r2_available(:, 2) - r2_targets(r2_i))); 
    end

    for j = 1: length(subj_indices)
        pattern = sprintf('vol_0*%d[^0-9]', subj_indices(j));
        loading_files(j) = filelist(cellfun(@(x) ~isempty(regexp(x, pattern)), filelist));
        r2_real(j) = r2_available(subj_indices(j), 2);
    end
    loading_files = fullfile(hyparams.folders.SUBJDIREFF, string(loading_files)');

    noneff_files = dir(hyparams.folders.SUBJDIRNONEFF);
    noneff_files = string({noneff_files(cellfun(@(x) contains(x, "vol_"), {noneff_files.name})).name});

    noneff_ind = randsample(length(noneff_files), length(loading_files));
    loading_files(:, 2) = fullfile(hyparams.folders.SUBJDIRNONEFF, noneff_files(noneff_ind))'; %this probably does not work!

    if any(p_noneffsubj > 0) % noneff effected subjects
        % noneff_indices = rand(length(loading_files), 1) < max(p_noneffsubj);
        n_additional_row = int32(length(r2_targets)*max(p_noneffsubj));
        noneff_ind = randsample(length(noneff_files), n_additional_row);
        loading_files(length(r2_targets) + 1 : length(r2_targets) + n_additional_row, :) = [ ...
            fullfile(hyparams.folders.SUBJDIRNONEFF, noneff_files(noneff_ind))',...
            loading_files(1 : n_additional_row, 2) ...
            ];
    end
end


function checketherogeneity(r2, nb_subj, nb_rep, subject_list, logger)
    r2_mean = []; r2_std = []; r2_all = []; r2_unique_mean = [];
    cumulative_r2 = [];
    for n = 1: nb_rep
        r2_targets = r2(1) + r2(2)*randn(nb_subj, 1);
        [loading_list, r2_n] = get_closest_samples( subject_list, ...
            r2_targets, file_list, opts.p_noneffsubj, hyparams);

        r2_unique_mean = [r2_unique_mean, length(unique(r2_n))];
        r2_mean = [r2_mean, mean(r2_n)];
        r2_std = [r2_std, std(r2_n)];
        r2_all = [r2_all, r2_n];
        cumulative_r2 = [cumulative_r2, r2_n];
    end
    logger.println("Target: " + r2(1) + " " + r2(1))
    logger.println("Subjects in the dataset: " + size(subject_list, 1) + "; nb of subj per iter: " + nb_subj)
    logger.println("Unique: mu: " + mean(r2_unique_mean) + "; std: " + std(r2_unique_mean) + "; max " + max(r2_unique_mean) + "; min: " + min(r2_unique_mean))
    logger.println("MEAN: mu: " + mean(r2_mean) + "; std: " + std(r2_mean) + "; max " + max(r2_mean) + "; min: " + min(r2_mean))
    logger.println("STDEV: mu: " + mean(r2_std) + "; std: " + std(r2_std) + "; max: " + max(r2_std) + "; min: " + min(r2_std))
    
    allR2Pickedup = unique(cumulative_r2);
    for r2_i = 1: length(allR2Pickedup)
        countR2Pickedup(r2_i) = sum(allR2Pickedup(r2_i) ==cumulative_r2);
    end
    logger.println("the 10 most occurent r2 effect subjects are repeated %s across 1000 sampling", ...
        string(maxk(countR2Pickedup, 10)).join());
end


function [hyparams, resfolder] = get_output_folders(hyparams, r2, datapath, opts, correction)
    if ~isempty(r2)
        hyparams.folders.r2_id = sprintf("r2_%s_%s", extractAfter(sprintf('%.3f', r2(1)) ,2), ...
            extractAfter(sprintf('%.3f', r2(2)), 2));
        if opts.vox_variability>0
            hyparams.folders.r2_id = hyparams.folders.r2_id + ...
                sprintf("wsv_%s_", extractAfter(sprintf('%.3f', opts.vox_variability), 2));
        end
    end

    % Extracting strings for identifying and create the results folder.
    expid = string( datetime( 'now', 'Format', 'yyMMdd_HHmmss'));
    grt_dst = ""; correction_char = "V"; p_noneffsubj=""; stdevnoneff=""; wsv="";
    if contains(datapath(1), "grt_dst"); grt_dst="grt_dst_"; end
    if contains(datapath(1), "wsv")
        pattern = 'wsv_(\d+)'; match = regexp(datapath(1), pattern, 'tokens', 'once');
        wsv=compose("wsv_%s_",match{1});
    end
    if correction == "clusterbased"; correction_char="K"; end
    if correction == "both"; correction_char="B"; end    
    if opts.p_noneffsubj > 0
        p_noneffsubj = sprintf( "pnes_%s_", extractAfter(sprintf('%.2f', opts.p_noneffsubj) ,2));
    end
    if ~isnan(opts.stdev_sample_fromnoneff)
        stdevnoneff = sprintf( "stdnoneff_%s_", extractAfter(sprintf('%.3f', opts.stdev_sample_fromnoneff) ,2));
    end
    
    % Folders for storing the results
    if strlength(opts.resfolder)==0
        resfolder = fullfile(hyparams.folders.RESULTS, grt_dst + ...
            hyparams.folders.r2_id + "_" + p_noneffsubj + stdevnoneff + wsv + correction_char + expid);
        hyparams.folders.LOG = fullfile(resfolder, 'log.txt');
        mkdir(resfolder);
    else
        resfolder = opts.resfolder;
        hyparams.folders.LOG = fullfile(resfolder, 'log.txt');
    end
end
