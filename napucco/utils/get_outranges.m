function [vx_range, vx_rangeII, subj_range] = get_outranges(INvx_range, ...
        INvx_rangeII, INsubj_range, hyparams)
    if isempty(INvx_range) % effected for voxelwise
        vx_range = 1;
    else; vx_range = INvx_range; end

    if isempty(INvx_rangeII) % NON-effected for voxelwise
        vx_rangeII = 500: 500: hyparams.second_cluster;
        vx_rangeII = unique([vx_rangeII, hyparams.second_cluster]);
    else; vx_rangeII = INvx_rangeII; end

    if isempty(INsubj_range)
        subj_range = 5:5:hyparams.nb_subj;
    else; subj_range = INsubj_range; end
end
