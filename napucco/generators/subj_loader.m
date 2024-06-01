function [p_values, R2s, eff_loader, noneff_loader] = subj_loader( ...
    hyparams, logger, subj_list, eff_loader, noneff_loader)
    % subj_list is a n_subj by 2 string array containing the files path to the
    % effected voxels and non-effected voxels. Loaded files must match with
    % the hyparams save.
    % If the number of voxels loaded is larger then the one specified in 
    % the parameters, only a first part of will be loaded.
    if length(subj_list) < hyparams.nb_subj; raiserrror; end
    nb_files = length(subj_list);
    addpath('./utils')
    addpath('./generators')

    % create subjects data and H0es and obtain a p-value matrix from them
    R2s.eff = cell(nb_files);
    R2s.noneff = cell(nb_files);
    p_values.eff = cell(nb_files, 1);
    p_values.noneff = cell(nb_files, 1);

    for subj_n = 1 : nb_files
        if (subj_list(subj_n, 1)~="none" && ~isfile(subj_list(subj_n, 1)))
            logger.clear()
            logger.println("Could not find %s ", subj_list(subj_n, 1), "tolog")
            RAISE ERROR
        end
        if (subj_list(subj_n, 2)~="none" && ~isfile(subj_list(subj_n, 2)))
            logger.clear()
            logger.println("Could not find %s ", subj_list(subj_n, 2), "tolog")
            RAISE ERROR
        end
        effpath = string(regexp(subj_list(subj_n, 1), '.*/(.*/[^/]+)$', 'tokens'));
        if isempty(effpath)
            effpath = "none";
        end
        noneffpath = string(regexp(subj_list(subj_n, 2), '.*/(.*/[^/]+)$', 'tokens'));
        if isempty(noneffpath)
            noneffpath = "none";
        end
        logger.print("loading effected: %s; non-eff: %s", ...
            effpath, ...
            noneffpath, "tolog")
        eff2load = subj_list(subj_n, 1);
        noneff2load = subj_list(subj_n, 2);
        if eff2load ~= "none"
            [data, eff_loader] = updatestack(eff_loader, eff2load);
            p_values.eff{subj_n} = data.plog_values(1:hyparams.numb_permutations + 1, 1: hyparams.nb_outfeat);
        end
        if noneff2load ~= "none"
            [data, noneff_loader] = updatestack(noneff_loader, noneff2load);
            p_values.noneff{subj_n} = data.plog_values(1:hyparams.numb_permutations + 1, 1: hyparams.nb_noneffvx);
        end
        logger.println("; loaded", "tolog");
    end

end

function [data, mstruct] = updatestack(mstruct, infile)
    if mstruct.init && ~isempty(mstruct.loaded_names) && any(infile == cat(1, mstruct.loaded_names{:}))
        ind2load = find(infile == cat(1, mstruct.loaded_names{:}));
        data = mstruct.data{ind2load};
        mstruct.stack_indices(mstruct.stack_indices==ind2load) = [];
        mstruct.stack_indices = [ind2load; mstruct.stack_indices];
    else
        data = load(infile); field = fieldnames(data);
        data = data.(field{1});
        if mstruct.init
            updateind = size(mstruct.data, 2) + 1;
            if updateind > mstruct.len
                updateind = mstruct.stack_indices(end); % otherwise retrive the last indices is the one less used
            end
            mstruct.stack_indices = [updateind; mstruct.stack_indices];
            mstruct.stack_indices = mstruct.stack_indices(1:min(size(mstruct.stack_indices,1), mstruct.len));
            mstruct.data{updateind} = data;
            mstruct.loaded_names{updateind} = infile;
        end
    end
end