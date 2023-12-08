function flattenpadding = encubator(a, b, n, m)
    %   ENCUBATOR  Wraps the effected voxels in a volume padded by
    %   non-effected voxels
    %   FLATTEN_WRAPPED_CUBE = ENCUBATOR(eff_data, padding_data, width_eff,
    %   padding_width) Given a effected data with width of the cube width_eff and 
    %   equally padding spaced padding_width towards every edge, it applies a
    %   cubic padding of the padding_data around the eff_data. It returns
    %   a flatten vector flattenpadding such that reshape(flattenpadding, 
    %   n, n, n) gives back the original cubic wrap.
    %

    flattenpadding = nan(length(a)+length(b), 1);
    ind_a = 1; ind_b = 1;
    ind_flattencube = 1;

    % A first m*(2*m+n)^2 volume of padding non-effected voxels is placed
    % before the part of the mid (in depth) part of the volume with
    % effected voxels
    flattenpadding(ind_flattencube:ind_flattencube+m*(2*m+n)^2-1) = b(ind_b:ind_b+m*(2*m+n)^2-1);
    ind_flattencube = ind_flattencube+m*(2*m+n)^2;
    ind_b = ind_b+m*(2*m+n)^2;

    for i=1:n % across the depth of the effected cube
        % A m*(2*m+n) surface is laid above (in height) the central part of the volume
        % (top of the sandwich)
        flattenpadding(ind_flattencube:ind_flattencube+m*(2*m+n)-1) = b(ind_b:ind_b+m*(2*m+n)-1);
        ind_flattencube = ind_flattencube+m*(2*m+n);
        ind_b = ind_b+m*(2*m+n);

        for j=1:n % across the height of the effected cube
            % m non effected voxels are placed ahead (in width - left side) in the width before
            % the effected ones
            flattenpadding(ind_flattencube:ind_flattencube+m-1) = b(ind_b:ind_b+m-1);
            ind_flattencube = ind_flattencube+m;
            ind_b = ind_b+m;
    
            % n effected voxels of one row along the width of the volume
            flattenpadding(ind_flattencube:ind_flattencube+n-1) = a(ind_a:ind_a+n-1);
            ind_flattencube = ind_flattencube+n;
            ind_a = ind_a+n;
    
            % m non effected voxels enclose the row (in width - right side)
            flattenpadding(ind_flattencube:ind_flattencube+m-1) = b(ind_b:ind_b+m-1);
            ind_flattencube = ind_flattencube+m;
            ind_b = ind_b+m;
    
        end

        % (bottom of the sandwich)
        flattenpadding(ind_flattencube:ind_flattencube+m*(2*m+n)-1) = b(ind_b:ind_b+m*(2*m+n)-1);
        ind_flattencube = ind_flattencube+m*(2*m+n);
        ind_b = ind_b+m*(2*m+n);
    end

    % Finally, a last m*(2*m+n)^2 volume of padding non-effected voxels is
    % placed at the end to wrap the volume
    flattenpadding(ind_flattencube:ind_flattencube+m*(2*m+n)^2-1) = b(ind_b:ind_b+m*(2*m+n)^2-1);
    ind_flattencube = ind_flattencube+m*(2*m+n)^2;
    ind_b = ind_b+m*(2*m+n)^2;
%     disp (ind_b + " " + length(b));
%     disp (ind_a + " " + length(a));
end