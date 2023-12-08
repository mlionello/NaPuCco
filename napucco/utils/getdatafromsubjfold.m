function [nb_files, filenames] = getdatafromsubjfold(folder)
    tmp = dir(folder);
    tmp = tmp(cellfun(@(x) ~isempty( ...
        regexp(x, 'subj_0*[0-9*].*\.mat$')), {tmp.name}));
    filenames = [];
    for i=  1:length(tmp)
        filenames = [filenames, string(tmp(i).name)];
    end
    nb_files = length(filenames);
end

