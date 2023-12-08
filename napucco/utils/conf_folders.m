function folders = conf_folders(nb_infeat, ...
    nb_timesteps, r2_target, grt_dist, vox_variability, nb_vx)
    if isempty(r2_target) || any(isnan(r2_target))
        r2_target = ["nan", "nan"];
    end
    folderID = compose("in%03d_t%04d",  nb_infeat,  nb_timesteps);
    exp_id = string( datetime( 'now', 'Format', 'yyMMdd_HHmmss'));
    r2_id = sprintf("r2_%s_%s", extractAfter(sprintf('%.3f',r2_target(1)) ,2), ...
    extractAfter(sprintf('%.3f',r2_target(2)), 2));
    vox_var_suffix = "";
    if vox_variability>0
        vox_var_suffix = sprintf("wsv_%s_", extractAfter(sprintf('%.3f',vox_variability) ,2));
        r2_id = r2_id + vox_var_suffix;
    end

    HOMEPATH = 'data';
    HOMEPATH = fullfile( HOMEPATH,  folderID);
    SUBJDIR = fullfile(HOMEPATH, 'subjects');

    if grt_dist == 0; SUBJDIREFF = fullfile(SUBJDIR, r2_id + "_" + exp_id);
    else; SUBJDIREFF = fullfile(SUBJDIR, 'grtdst_'  + vox_var_suffix + exp_id); end
    SUBJDIRNULL = fullfile(SUBJDIR,  "nulldst_" + exp_id);
    RESULTS = fullfile(HOMEPATH, 'results');
    LOGS = fullfile(HOMEPATH, 'logs');
    FIGURES = fullfile(HOMEPATH, 'figures');

    if ~ isfolder(HOMEPATH); mkdir(HOMEPATH); end
    if ~ isfolder(SUBJDIR); mkdir(SUBJDIR); end
    if nb_vx(2) == 0; SUBJDIRNULL="none"; else; mkdir(SUBJDIRNULL); end
    if nb_vx(1) == 0; SUBJDIREFF="none"; else; mkdir(SUBJDIREFF); end
    if ~ isfolder(FIGURES); mkdir(FIGURES); end
    if ~ isfolder(LOGS); mkdir(LOGS); end
    if ~ isfolder(SUBJDIR); mkdir( SUBJDIR); end
    if ~ isfolder(RESULTS); mkdir( RESULTS); end

    folders = struct("folderID", folderID, 'exp_id', exp_id, ...
        'r2_id', r2_id, ...
        "HOMEPATH", HOMEPATH, ...
        "SUBJDIR", SUBJDIR, ...
        "SUBJDIREFF", SUBJDIREFF, ...
        "SUBJDIRNULL", SUBJDIRNULL, ...
        "LOGS", LOGS, "RESULTS", RESULTS, "FIGURES", FIGURES);
end
