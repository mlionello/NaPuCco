function superman = nonpc(p_values, hyparams, subj_range, n_additional_row, p_noneff_range, nb_cpus)

    nb_subj_range = length(subj_range);
    nb_noneff_subj_range = length(p_noneff_range); % p_noneff must at least be [0]

    tmpeff = nan(nb_subj_range*nb_noneff_subj_range, ...
        hyparams.numb_permutations + 1, size(p_values.eff{1}, 2));
    tmpnoneff = nan(nb_subj_range*nb_noneff_subj_range, ...
        hyparams.numb_permutations + 1, size(p_values.noneff{1}, 2));

    tot_eff_subjs = numel(p_values.eff) - n_additional_row;

    cr = ""; tic;
    nb_iter = nb_subj_range*nb_noneff_subj_range;
    parfor_progress(nb_iter);
    for k = 1 : nb_iter
        if ~isempty(getCurrentTask())
            parfor_progress;
        else
            msg = compose(" subj range %d / %d ", k, nb_subj_range*nb_noneff_subj_range);
            fprintf("" + cr + msg)
            cr = repmat('\b', 1, strlength(msg));
        end
        [j, i] = ind2sub([nb_subj_range, nb_noneff_subj_range], k);
        nb_subjs_in = subj_range(j);
        nb_noneff_in = int32(p_noneff_range(i)*nb_subjs_in);

        onset_init = tot_eff_subjs - nb_subjs_in;
        final_onset = onset_init + nb_noneff_in + 1;
        final_offset = final_onset + subj_range(j) - 1;

        indices = final_onset: final_offset;

        tmpeff(k, :, :) = combining_funct(p_values.eff( indices));

        tmpnoneff(k, :, :) = combining_funct(p_values.noneff( indices ));
    end
    parfor_progress(0); toc;

    superman.eff = tmpeff;
    superman.noneff = tmpnoneff;

end


% a = 1:14;
% disp(a)
% p_noneff = [0, 0.2, 0.5];
% tot_eff_subjs = 9;
% subj_range = [3:2:9];
% for k = 1 : numel(p_noneff)*numel(subj_range)
%         [j, i] = ind2sub([numel(subj_range), numel(p_noneff)], k);
%         nb_subjs_in = subj_range(j);
%         offset_p_noneff = int32(p_noneff(i)*nb_subjs_in);
%         onset_init = tot_eff_subjs - nb_subjs_in;
%         final_onset = onset_init + offset_p_noneff + 1;
%         final_offset = final_onset + subj_range(j) - 1;
% 
%         fprintf("%d ", a(final_onset:final_offset))
%         fprintf("\n")
% end