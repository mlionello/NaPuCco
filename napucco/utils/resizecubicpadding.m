function [a_indices, b_indices, n, m] = resizecubicpadding(a, b, logger)
    arguments
        a;
        b;
        logger = nan;
    end
    n = cubicrt(a);
    a_indices = 1:int32(n^3);
    m = get_paddingwidth(n, b);
    if (2*int32(m)+n)^3==(b+length(a_indices)); m=int32(m); end
    m = floor(m);
    b_indices = 1: (n+2*m)^3-n^3;
%     if isnan(logger)
%         fprintf('first array is %d long (cubrt %.2f), the second one is %d long', ...
%         a, n, b)
%         fprintf('maximum padding available is the floor of %.2f',  m)  
%         fprintf(['final size of effected voxels: %d (cubrt %d), ' ...
%         'non-effected voxels: %d, (padding %d)'], length(a_indices), n, ...
%         length(b_indices), m)
%     else
        logger.println('first array is %d long (cubrt %.2f), the second one is %d long', ...
            a, n, b)
        logger.println('maximum padding available is the floor of %.2f',  m)  
        logger.println(['final size of effected voxels: %d (cubrt %d), ' ...
        'non-effected voxels: %d, (padding %d)'], length(a_indices), n, ...
        length(b_indices), m)

%     end
end
