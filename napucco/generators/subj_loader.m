function [p_values, R2s, eff_loader, null_loader] = subj_loader( ...
    hyparams, logger, subj_list, eff_loader, null_loader)
    % subj_list is a n_subj by 2 string array containing the files path to the
    % effected voxels and non-effected voxels. Loaded files must match with
    % the hyparams save.
    % If the number of voxels loaded is larger then the one specified in 
    % the parameters, only a first part of will be loaded.
    if length(subj_list) ~= hyparams.nb_subj; raiserrror; end
    addpath('./utils')
    addpath('./generators')

    % create subjects data and H0es and obtain a p-value matrix from them
    R2s.eff = cell(hyparams.nb_subj);
    R2s.null = cell(hyparams.nb_subj);
    p_values.eff = cell(hyparams.nb_subj, 1);
    p_values.null = cell(hyparams.nb_subj, 1);

    for subj_n = 1 : hyparams.nb_subj
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
        nullpath = string(regexp(subj_list(subj_n, 2), '.*/(.*/[^/]+)$', 'tokens'));
        if isempty(nullpath)
            nullpath = "none";
        end
        logger.print("loading effected: %s; non-eff: %s", ...
            effpath, ...
            nullpath, "tolog")
        eff2load = subj_list(subj_n, 1);
        null2load = subj_list(subj_n, 2);
        if eff2load ~= "none"
            [data, eff_loader] = updatestack(eff_loader, eff2load);
            p_values.eff{subj_n} = data.p_values(1:hyparams.numb_permutations + 1, 1: hyparams.nb_outfeat);
            R2s.eff{subj_n} = data.R2s(1: hyparams.nb_outfeat);
        end
        if null2load ~= "none"
            [data, null_loader] = updatestack(null_loader, null2load);
            p_values.null{subj_n} = data.p_values(1:hyparams.numb_permutations + 1, 1: hyparams.second_cluster);
            R2s.null{subj_n} = data.R2s(:, 1: hyparams.second_cluster);
        end
        logger.println("; loaded", "tolog");
    end
    p_values.eff = permute(cat(3, p_values.eff{:}), [3 1 2]);
    p_values.null = permute(cat(3, p_values.null{:}), [3 1 2]);
    R2s.eff = permute(cat(2, R2s.eff{:}), [2 1]);
    R2s.null = permute(cat(2, R2s.null{:}), [2 1]);
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