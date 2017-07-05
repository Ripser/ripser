function [Is] = ripserPC( X, coeff, maxdim, thresh )
    %:param X: N X d point cloud with N points in d dimensions
    %:param coeff: Field coefficient with which to run TDA
    %:param maxdim: Maximum dimension of homology
    %:param thresh: Threshold up to which to add edges
    D = pdist(X);
    if nargin < 4
        thresh = max(D(:))*2;
    end
    d = single(D(:));
    Is = ripser(d, coeff, maxdim, thresh);
end

