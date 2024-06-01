function [nvox, clbc_n_range, mvox, clbc_m_range] = check_clusterbased_padding( ...
            nb_outfeat, second_cluster, clbc_n_range, clbc_m_range, correction_char, opts, logger)
    % Checking for voxels qunatity for clusterbased-correction
    eff_cube_edge = cubicrt(nb_outfeat);
    nvox = int32(eff_cube_edge^3); mvox = nb_outfeat;
    % If there are have less voxels than the ones required in the range, then reduce the
    % range
    if eff_cube_edge < clbc_n_range(end)
        clbc_n_range = [1: eff_cube_edge];
        nvox = int32(eff_cube_edge^3);
        eff_cube_edge = cubicrt(nvox);
        logger.println("The amount of voxels availables %d are not enough for the " + ...
            "maximum cubic range. Only %d will be used to match a cube of size %d ^3", ...
            nb_outfeat, nvox, clbc_n_range(end));
    end
    % If there are more voxels than the ones required in the range, then reduce the
    % number of given voxels to match the maximum number of voxels in the range
    if eff_cube_edge > clbc_n_range(end)
        nvox = int32(clbc_n_range(end)^3);
        eff_cube_edge = cubicrt(nvox);
        logger.println("The amount of voxels availables %d are exceeding the ones " + ...
            "needed by the maximum cubic range. Only %d will be used to match a cube of size %d \^3", ...
            nb_outfeat, nvox, clbc_n_range(end));
    end
    null_padding_width = int32(get_paddingwidth( eff_cube_edge, mvox));
    % If there are have less null voxels than the ones required in the range, then reduce the
    % number of maximum voxels in the range to match the null voxels
    if opts.fpr && mvox < (2*clbc_m_range(end) + eff_cube_edge)^3-eff_cube_edge^3
        clbc_m_range = [1: null_padding_width];
        mvox = (2*clbc_m_range(end) + eff_cube_edge)^3 - eff_cube_edge^3;
        logger.println("The amount of null voxels availables %d are not enough for the " + ...
            "maximum PADDING cubic range. Only %d will be used to match a PADDING of size %d", ...
            second_cluster, mvox, clbc_m_range(end));
    end
    % If there are have more null voxels than the ones required in the range, then reduce the
    % number of voxels  to match the maximum number of voxels in the range
    if opts.fpr &&  mvox > (2*clbc_m_range(end) + eff_cube_edge)^3-eff_cube_edge^3
        mvox = (2*clbc_m_range(end) + eff_cube_edge)^3 - eff_cube_edge^3;
        logger.println("The amount of null voxels availables %d are exceeding the ones " + ...
            "needed by the PADDING maximum cubic range. Only %d will be used to match a PADDING of size %d", ...
            second_cluster, mvox, clbc_m_range(end));
    end
end