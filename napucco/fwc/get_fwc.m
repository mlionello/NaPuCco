function [p_fwc, H0_maxima, p_uncorrected] = get_fwc(A)
    % GET_FWC compute familywise corrected values.
    %   [p_fwc, H0_maxima] = GET_FWC(A) computes the familiwise corrected
    %   values for the first row of A with respect to the other rows.
    %   It returns the fwc-values for the first row and a vector of maxima
    %   collected from the A sub-matrix (2:end, :) for each row.

    % A input matrix can be u_values: [J+1, L]; p_values: [N, J+1, L];
    % [J+1, L]      ->  [L]
    % [N, J+1, L]   ->  [L; N]
    % FWC:

    %prepare dims match whether input matrix has either ndims 2 (u-values)
    %or 3 (p-values). A becomes to match : [nb_perm + 1, nb_subj, nb_voxels]
    p_uncorrected = nan;
    if isempty(A)
        p_fwc = nan; H0_maxima = nan; p_uncorrected = nan; return
    end

    if ndims(A) == 3
        A = permute(A, [2 1 3]);
    elseif ndims(A) == 2
        A = reshape(A, size(A, 1), 1, size(A, 2));
    end
    %VoxelWise correction: for each permuation (first dim) extract the maximum 
    %value across the voxels (last dim)
    H0_maxima = max(A(2:end, :, :), [], 3); % MAXIMUM IF FISHER

    %Family-wise correction: compute the family-wise corrected values for
    %each voxel in the first row with respect the maxima extracted in
    %the previous step

    repH0 = repmat(H0_maxima, [1, 1, size(A, ndims(A))]);
    mixed = cat(1, A(1, :, :), repH0); % same of doing [~, I] = sort(mixed)?
    fwc = tiedrank( -mixed) /size( mixed, 1); % MINUS if fisher
    p_fwc = reshape(fwc(1, :, :), size(fwc, 2), size(fwc, 3));
    p_fwc = p_fwc';
    % p_fwc = permute(fwc(1, :, :), [3, 2 ,1]);

    %p_uncorrected = tiedrank( -A) /size( A, 1);
    %p_uncorrected = squeeze(p_uncorrected(1, :, :));
    %p_uncorrected = p_uncorrected';

%     tmp = squeeze(A);
%     hist(max(tmp(2:end,:),[],2))
%     tmp2 = sort(tmp(1,:));
end
