function [rates, tot_vx_thr, tot_vx] = get_rates(voxels, pies_thr)
    rates = containers.Map();
    for p_i = 1: length(pies_thr)
        tot_vx_thr = squeeze(sum(voxels<pies_thr(p_i)));
        tot_vx = length(voxels)-squeeze(sum(isnan(voxels)));
        rates(string(pies_thr(p_i))) = tot_vx_thr./tot_vx;
    end
end
