function fwc_subjextraction(homepath, nb_subj, nb_rep, opts)
arguments
    homepath;
    nb_subj = 20;
    nb_rep = 100;
    opts.power = 1;
    opts.fpr = 1;
end
    addpath('./utils')
    addpath('./fwc')
    addpath('./generators')
    
    %main(nb_subj, nb_infeat, nb_outfeat, nb_timesteps, rsquared_target, correction, opts)
    homepath = "data/" + homepath;
        
    for i = 1: nb_rep
            [fwc, ~] = voxel_cluster_correction(...
                "voxelwise", ...
                nb_subj, ...
                'saving', 0, ...
	            'save_res', 0, ...
                'homepath', homepath);
    
	    if ~isdir(fullfile(homepath, "results", "grouping"))
		    mkdir(fullfile(homepath, "results", "grouping")); end
        if opts.power
            fwctmp = fwc.results_fwc_uvalues.power('0.05');
            save(fullfile(homepath, "results", "grouping", "power05_"+compose("%03d",i)+".mat"), "fwctmp");
            fwctmp = fwc.results_fwc_uvalues.power('0.001');
            save(fullfile(homepath, "results", "grouping", "power001_"+compose("%03d",i)+".mat"), "fwctmp");
            fwctmp = fwc.results_fwc_uvalues.power('0.0001');
            save(fullfile(homepath, "results", "grouping", "power0001_"+compose("%03d",i)+".mat"), "fwctmp");
        end
        if opts.fpr
            fwctmp = fwc.results_fwc_uvalues.fpr('0.05');
            save(fullfile(homepath, "results", "grouping", "fpr05_"+compose("%03d",i)+".mat"), "fwctmp");
            fwctmp = fwc.results_fwc_uvalues.fpr('0.001');
            save(fullfile(homepath, "results", "grouping", "fpr001_"+compose("%03d",i)+".mat"), "fwctmp");
            fwctmp = fwc.results_fwc_uvalues.fpr('0.0001');
            save(fullfile(homepath, "results", "grouping", "fpr0001_"+compose("%03d",i)+".mat"), "fwctmp");
        end
    end            
    
end
