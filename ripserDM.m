function [Is] = ripserDM( D, coeff, maxdim, thresh, wrap )
    %:param D: (N X N distance matrix
    %:param coeff: Field coefficient with which to run TDA
    %:param maxdim: Maximum dimension of homology
    %:param thresh: Threshold up to which to add edges
    %:param wrap: Whether to wrap around ripser with a system call
    % (for people having trouble with mex who can compile locally)
    if nargin < 4
        thresh = max(D(:))*2;
    end
    if nargin < 5
        wrap = 0;
    end
    
    if wrap
        dlmwrite('D.txt', D);
        cmd = sprintf('./ripser-coeff --format distance --dim %i --modulus %i --threshold %g D.txt', maxdim, coeff, thresh);
        disp(cmd)
        system(cmd); 
    else
        %pdist surrogate
        N = size(D, 1);
        [I, J] = meshgrid(1:N, 1:N);
        d = D(I < J);
        d = single(d(:));
        Is = ripser(d, coeff, maxdim, thresh);
    end
    
    %Delete entries with numerically insignificant persistence
    for ii = 1:length(Is)
        I = Is{ii};
        if numel(I) == 0
            continue;
        end
        P = I(:, 2) - I(:, 1);
        Is{ii} = I(P > eps, :);
    end
end

