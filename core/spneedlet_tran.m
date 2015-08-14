function beta = spneedlet_tran(alm, l_max, B)
%SPNEEDLET_TRAN   Fast spherical needlet transform. It computes the 
%needlet coefficients based on the spherical harmonic coefficients alm.
%
%   spneedlet_tran(alm, l_max, B) 
%
% Inputs:
%   alm - the spherical harmonic coefficients, (l_max+1)-by-(2*l_max+1)
%   matrix, alm(l+1, m+l_max+1) is the spherical harmonic coefficient with
%   subscripts l and m
%   l_max - the maximal l
%   B - the parameter
%
% Outputs:
%   beta - the needlet coefficients, (j_max+1)-by-1 cell, the (j+1)-th cell
%   is the needlet coefficients with subscript j
%
% Author: Minjie Fan, 2015

t = cputime;

bw = l_max+1;

j_max = get_j_max(B, l_max);

beta = cell(j_max+1, 1);

b_vector = get_b_vector(B, j_max, l_max);

blmj = zeros(size(alm));
    
for j = 0:j_max
    disp(['j = ', num2str(j), ' starts...']);
    Nside = get_Nside(B, j);
    disp(['Nside = ', num2str(Nside)])
    
    Npix = 12*Nside^2;
    
    legendre_norm = get_legendre_norm(l_max, Nside);

    l_st = ceil(B^(j-1));
    l_en = min(floor(B^(j+1)), l_max);
    
    for l = l_st:l_en
        blmj(l+1, bw:l+bw) = alm(l+1, bw:l+bw)*b_vector(j+1, l);
    end
    
    beta{j+1} = inv_spharmonic_tran(blmj, Nside, legendre_norm, bw, l_st, l_en)*...
        sqrt(4*pi/Npix);
   
end

e = cputime-t;
disp(['Elapsed CPU time: ', num2str(e)])

end
    
    
    