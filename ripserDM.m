function [Is] = ripserDM( D, coeff, maxdim, thresh )
    if nargin < 4
        thresh = max(D(:))*2;
    end
    d = squareform(D);
    d = single(d(:));
    Is = ripser(d, coeff, maxdim, thresh);
end

