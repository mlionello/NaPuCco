classdef Napucco
    % NaPuCco - A MATLAB class for performing group-level inference of 
    % fMRI data through non-parametric combination
    %
    %   This class provides functionality for generating volumes, running
    %   experiments, and computing results.
    %
    %   Usage:
    %       napuccoObj = Napucco(nb_infeat, nb_timepoints)
    %       napuccoObj.fetch_volumes(folders)
    %       napuccoObj.generate_volumes(nb_volumes, opts)
    %       napuccoObj.run_experiments(correction, r2, nb_subj, nb_rep, opts)
    %       Napucco.compute_result()
    %
    %  Examples:
    %       gli_npc = Napucco(nb_infeat, nb_timepoints)
    %       gli_npc = gli_npc.generate_volumes(20, 'nb_effvx', 27, 'nb_noneffvx', 1000, 'r2_target', [0.07, 0.02]);
    %       gli_npc = gli_npc.gli_npc.run_experiments('both', [])
    %       gli_npc.compute_result()
    %
    
    properties
        nb_timepoints
        nb_infeat
        folders
        setup
        results_path
    end

    methods
        % Constructor
        function obj = Napucco(nb_infeat, nb_timepoints)
            % Napucco - Constructor for the Napucco class.
            %
            %   obj = Napucco(nb_infeat, nb_timepoints)
            %
            %   Parameters:
            %       nb_infeat (int): Number of regressors.
            %       nb_timepoints (int): Number of timepoints.

            obj.nb_infeat = nb_infeat;
            obj.nb_timepoints = nb_timepoints;
        end

        % Method to fetch volumes if separatedly computed
        function obj = fetch_volumes(obj, folders)
            % fetch_volumes - Fetches volumes from specified folders.
            %
            %   obj = fetch_volumes(obj, folders)
            %
            %   Parameters:
            %       folders (struct): Structure containing folder paths.
            %           - eff_path (string): Path for effective volumes.
            %           - noneff_path (string): Path for non-effective volumes.

            arguments
                obj;
                folders.eff_path = [];
                folders.noneff_path = [];
            end

            if ~isempty(folders.eff_path)
                obj.folders.eff = folders.eff_path;
            end
            if ~isempty(folders.noneff_path)
                obj.folders.noneff = folders.noneff_path;
            end
        end

        % Method to generate volumes
        function obj = generate_volumes(obj, nb_volumes, opts)
            % generate_volumes - Generates volumes based on specified options.
            %
            %   obj = generate_volumes(obj, nb_volumes, opts)
            %
            %   Parameters:
            %       nb_volumes (int): Number of volumes to generate.
            %       opts (struct): Options for volume generation.
            %           - nb_effvx (int): Number of effective voxels.
            %           - nb_noneffvx (int): Number of non-effective voxels.
            %           - r2_target: Target R2 value.
            %           - grtdst (int): Greater distance.
            %           - path (string): Path for resuming generation.

            arguments
                obj;
                nb_volumes int32;
                opts.nb_effvx int32 = 0;
                opts.nb_noneffvx int32 = 0;
                opts.r2_target = [];
                opts.grtdst int32 = 0;
                opts.path string = [];
            end

            if ~isempty(opts.path)
                resume_generation(opts.path);
            else
                [hyparams, ~] = vc_correction('none', ...
                    nb_volumes, obj.nb_infeat, opts.nb_effvx, obj.nb_timepoints, ...
                    opts.r2_target, 'sec_cluster_voxels', opts.nb_noneffvx, ...
                    'greater_dist', opts.grtdst);
            end

            if opts.nb_noneffvx > 0
                obj.folders.noneff = hyparams.folders.SUBJDIRNULL;
            end
            if opts.nb_effvx > 0
                obj.folders.eff = hyparams.folders.SUBJDIREFF;
            end
            obj.setup = false;
        end

        % Method to run experiments
        function obj = run_experiments(obj, correction, r2, nb_subj, nb_rep, opts)
            % run_experiments - Simulates experiments.
            %
            %   obj = run_experiments(obj, correction, r2, nb_subj, nb_rep, opts)
            %
            %   Parameters:
            %       correction (string): Type of correction ('voxelwise', 'clusterbased', 'both').
            %       r2: R2 value.
            %       nb_subj (int): Number of subjects.
            %       nb_rep (int): Number of repetitions.
            %       opts (struct): Additional options.
            %           - fpr (bool): False positive rate.
            %           - power (bool): Power.
            %           - resfolder: Result folder.

            arguments
                obj;
                correction string {mustBeMember(correction, {'voxelwise', 'clusterbased', 'both'})};
                r2;
                nb_subj int32 = 80;
                nb_rep int32 = 100;
                opts.fpr int32 = 1;
                opts.power int32 = 1;
                opts.resfolder = nan;
            end

            if isempty(obj.folders)
                error('volumes feeding first');
            end

            obj.results_path = subj_subsampling( ...
                correction, [obj.folders.eff, obj.folders.noneff], nb_subj, ...
                nb_rep, r2, 'fpr', opts.fpr, 'power', opts.power, ...
                'resfolder', opts.resfolder);
        end
    end

    methods(Static)
        % Static method to compute results
        function compute_result(obj)
            % compute_result - Computes and analyzes simulation results.
            %
            %   Napucco.compute_result()
            
            if isempty(obj.results_path)
                error('simulate some (>100) experiments first!');
            end

            for corr_method = ["voxelwise_correction", "clusterbased_correction"]
                res_subfolder = fullfile(obj.results_path, corr_method);
                if exist(res_subfolder, 'dir')
                    compute_mean_res(res_subfolder);
                end
            end
        end
    end
end
