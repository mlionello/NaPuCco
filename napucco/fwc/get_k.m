function K = get_k(mask, v)
arguments
    mask;
    v=1;
end

    N = int32(length(mask).^(1/3));
    visited = ~mask;
    K = 0;
    if exist('bwconncomp')
        bb = bwconncomp(reshape(mask, [N,N,N]), 18);
        numPixels = cellfun(@numel,bb.PixelIdxList);
        [K,~] = max(numPixels);
        if isempty(K)
            K=0;
        end
    else
        if v
            for i = 1:length(mask)
                if ~visited(i)
                    [K(end+1), visited] = browse_neigh_flat(visited, i, N);
                end
            end
        else
            visited = reshape(visited, N, N, N);
            for z = 1: N
                for y = 1: N
                    for x = 1: N
                        if ~visited(x, y, z)
                            [K(end+1), visited] = browse_neigh(visited, [x, y, z]);
                        end
                    end
                end
            end
        end
    end
end

function [k, visited] = browse_neigh(visited, coords)
    if any(coords>(size(visited)), 'all') || ...
        any(coords<1, 'all') || ...
        visited(coords(1), coords(2), coords(3)) %visited(coords) might be omitted
        k = 0;
        return
    end
    k(1) = 1;
    visited(coords(1), coords(2), coords(3)) = 1;
    [k(2), visited] = browse_neigh(visited, coords-int32([1,0,0]) );
    [k(3), visited] = browse_neigh(visited, coords+int32([1,0,0]) );
    [k(4), visited] = browse_neigh(visited, coords-int32([0,1,0]) );
    [k(5), visited] = browse_neigh(visited, coords+int32([0,1,0]) );
    [k(6), visited] = browse_neigh(visited, coords-int32([0,0,1]) );
    [k(7), visited] = browse_neigh(visited, coords+int32([0,0,1]) );
    k = sum(k);
end

function [k, visited] = browse_neigh_flat(visited, i, N)
    if i>numel(visited) || i<1 || visited(i)
        k = 0; 
        return
    end
    k(1) = 1;
    visited(i) = 1;
    if mod(i-1, N^2)+N+1<=N^2
        [k(2), visited] = browse_neigh_flat(visited, i+N, N);
    end
    if mod(i-1, N^2)-N+1>=1
        [k(3), visited] = browse_neigh_flat(visited, i-N, N);
    end
    if mod(i-1, N^3)+N^2+1<=N^3
        [k(4), visited] = browse_neigh_flat(visited, i+N^2, N);
    end
    if mod(i-1, N^3)-N^2+1>=1
        [k(5), visited] = browse_neigh_flat(visited, i-N^2, N);
    end
    if mod(i-1, N)+1+1<=N
        [k(6), visited] = browse_neigh_flat(visited, i+1, N);
    end
    if mod(i-1 ,N)-1+1>=1
        [k(7), visited] = browse_neigh_flat(visited, i-1, N);
    end
    k = sum(k);
end
