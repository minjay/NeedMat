function beta = spneedlet_tran(alm, l_max, B)
%SPNEEDLET_TRAN Needlet transform.
%
%   spneedlet_tran( coef, l_max, B ) 
%
% Inputs:
%   coef - the spherical harmonic coefficients. coef[i, j] gives a[l,
%   m]=a[i-1, j-l_max-1]
%   l_max - the maximal value of l
%   B - the parameter
% Outputs:
%   beta - the needlet coefficients

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
    
    
    