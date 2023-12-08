function sorted_indices = cubic_radial_sort(nb_voxels)
%CUBIC_RADIAL_SORT It returns the indices to apply to a descend sorted
% 1D-array and get a radial distribution if reshaped into a cube
% e.g.    [R2s, outdata_indices] = sort(R2s, "descend");
%         R2s = R2s(sorted_indices);
%         radialsortedcube = reshape(R2s, cubicrt(nb_voxels), cubicrt(nb_voxels), cubicrt(nb_voxels))
%         pv_param = pv_param(outdata_indices(sorted_indices));
        cube_size = double(cubicrt(nb_voxels));
        [x, y, z] = meshgrid(1:cube_size, 1:cube_size, 1:cube_size);
        center_point = ceil((cube_size + 1) / 2);
        radial_distance = sqrt((x - center_point).^2 + (y - center_point).^2 + (z - center_point).^2);
        flatten_indices = reshape(radial_distance, 1, []);
        [~, indices] = sort(flatten_indices);
        [~, sorted_indices] = sort(indices);
end

