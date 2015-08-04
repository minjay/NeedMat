function Spharmonic_coef = spharmonic_tran( fun_sample_file, bw )
%SPHARMONIC_TRAN Spherical harmonic transform using the semi-naive
%DLT. Precomputes all necessary associated Legendre functions prior to
%transforming.
%
%   spharmonic_tran( fun_sample_file, bw ) is the spherical
%   harmonic transform of the function samples stored in the file with name
%   fun_sample_file. The samples in the file fun_sample_file are arranged
%   in interleaved real/imaginary format. bw denotes the bandwidth of
%   the bandlimited function. SPHARMONIC_TRAN returns the coefficients in
%   the bw by 2bw-1 matrix Spharmonic_coef. The entry in the ith row and the
%   jth column represents the coefficient corresponding to l=i-1, m=j-B.
%
%   See the packages SpharmonicKit and S2kit for details.
%
%   Written by Minjie Fan

filename = ['coef', fun_sample_file];
system(['./test_s2_semi_memo_for ', fun_sample_file, ' ', filename, ' ',...
    num2str(bw),' 0']);
tmp = textread( filename );
index = 0;
Spharmonic_coef = zeros(bw, 2*bw-1);
for m = 0:bw-1
    for l = m:bw-1
        index = index+1;
        Spharmonic_coef(l+1,m+bw) = tmp( 2*index-1 )+tmp( 2*index )*1i;
    end
end
for m = 1-bw:-1
    for l = abs(m):bw-1
        index = index+1;
        Spharmonic_coef(l+1,m+bw) = tmp( 2*index-1 )+tmp( 2*index )*1i;
    end
end
system(['rm ', filename]);

end