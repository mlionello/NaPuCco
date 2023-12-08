function [hyparams,  dist_settings] = ...
        setparameters(nb_outfeat, opts, nb_subj, nb_timesteps, nb_infeat, ...
        rsquared_target, folders)

    subj_mean_range = opts.subj_mean_range;
    subj_std_range = opts.subj_std_range;
    threshold = opts.threshold;
    numb_permutations = opts.numb_permutations;
    vox_std_range = opts.vox_std_range;

    shuffler = nan(numb_permutations, nb_timesteps);
    for perm = 1: numb_permutations
        shuffler(perm, : ) = randperm( nb_timesteps);
    end
    shufflernull = nan(numb_permutations, nb_timesteps);
    for perm = 1: numb_permutations
        shufflernull(perm, : ) = randperm( nb_timesteps);
    end
    
    hyparams = struct("nb_subj", nb_subj, "nb_timesteps", nb_timesteps, ...
        "nb_infeat", nb_infeat, "nb_outfeat", nb_outfeat, ...
        "numb_permutations", numb_permutations, "second_cluster", opts.sec_cluster_voxels, ...
        "shuffler", shuffler, "shufflernull", shufflernull, "folders", folders);
    dist_settings = struct("rsquared_target", rsquared_target, "threshold", threshold, ...
        "subj_mean_range", subj_mean_range, "subj_std_range", subj_std_range, ...
        "vox_std_range", vox_std_range, "vox_variability", opts.vox_variability);
    
    disp(array2table(struct2cell(hyparams)','VariableNames', fieldnames(hyparams)))
    disp(array2table(struct2cell(dist_settings)','VariableNames', fieldnames(dist_settings)))
    
end
