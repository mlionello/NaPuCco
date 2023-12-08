function superman = nonpc(p_values, hyparams, subj_range)
    tmpeff = nan(length(subj_range), hyparams.numb_permutations + 1, size(p_values.eff, 3));
    tmpnull = nan(length(subj_range), hyparams.numb_permutations + 1, size(p_values.null, 3));

    %p_values = p_values(randperm(size(p_values, 1)), :, :);
    for j = 1: length(subj_range)
        tmpeff(j, :, :) = combining_funct(p_values.eff( 1: subj_range(j), :, :));
        tmpnull(j, :, :) = combining_funct(p_values.null( 1: subj_range(j), :, :));
    end

    superman.eff = tmpeff;
    superman.null = tmpnull;

end

