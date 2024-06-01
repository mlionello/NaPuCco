function resume_generation(path)
arguments
    path;
end
% it reasumes the generation with the parameters obtained from the given
% path
    addpath('./utils')
    addpath('./generators')
    addpath('./fwc')
        
    path = fullfile( path);
    [dist_settings, hyparams] = parse_settings(path);

    generate_volumes( ...
        hyparams.nb_subj, hyparams.nb_infeat, hyparams.nb_timesteps, ...
        hyparams.nb_outfeat, dist_settings.rsquared_target, ...
        'prev_settings', struct('hyparams', hyparams, ...
            'dist_settings', dist_settings));
end

function [dist_settings, hyparams] =  ...
        parse_settings(path)
    dist_settings = load(fullfile(path, "dist_settings.mat"));
    dist_settings = dist_settings.dist_settings;
    hyparams = load(fullfile(path, "hyparams.mat"));
    hyparams = hyparams.hyparams;
end