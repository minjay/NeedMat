function map = inv_spharmonic_tran_naive(alm, theta, phi, l_max)
%INV_SPHARMONIC_TRAN_NAIVE   A naive implementation of the inverse 
%spherical harmonic transform.
%
%   map = inv_spharmonic_tran_naive(alm, theta, phi, l_max)
%
% Inputs:
%   alm - the spherical harmonic coefficients, (l_max+1)-by-(2*l_max+1)
%   matrix, alm(l+1, m+l_max+1) is the spherical harmonic coefficient with
%   subscripts l and m
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%   l_max - the maximal l
%
% Outputs:
%   map - the reconstructed map, N-by-1 vector
%
% Author: Minjie Fan, 2015

n = length(theta);
map = zeros(n, 1);
for i = 1:n
    if mod(i, 1000)==0
        i
    end
    for l = 0:l_max
        for m = 0:l
            if m==0
                map(i) = map(i)+alm(l+1, m+l_max+1)*spharmonic_eval(l, m, theta(i), phi(i));
            else
                map(i) = map(i)+2*real(alm(l+1, m+l_max+1)*spharmonic_eval(l, m, theta(i), phi(i)));
            end
        end
    end
end

end
