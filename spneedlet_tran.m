function beta = spneedlet_tran( fun_sample_file, l_max, B )
%SPNEEDLET_TRAN Needlet transform.
%
%   spneedlet_tran( fun_sample_file, l_max, B ) does a needlet transform of
%   the function whose samples are given in the file with name
%   fun_sample_file. l_max+1 is the bandwidth of the function. B (B>1)
%   represents the parameter in the needlet transform. The return value
%   beta gives the spherical needlet coefficients.
%
%   This function depends on spharmonic_tran, inv_spharmonic_tran and pix2ang in the package
%   MEALPix.
%
%   Written by Minjie Fan

t = cputime;

%compute j_max from the inequality B^{j-1}<=l_max
j_max = fix( log(l_max)/log(B)+1 );
if abs(j_max-(log(l_max)/log(B)+1))<eps
    j_max = j_max-1;
end

%do the spherical harmonic transform
bw = l_max+1;
coef = spharmonic_tran( fun_sample_file, bw );

%do the needlet transform, in which we:
%(1)    evaluate the window function b before doing the transform
%(2)    compute cubature points and cubature weights for j=0,...,j_max by using
%the MEALPix package

%beta records the spherical needlet coefficients
beta = cell(j_max+1,1);

%b_vector records the evaluations of the window function b
b_vector = zeros(j_max+1,l_max);
for j = 0:j_max
    for l = 1:l_max
        b_vector(j+1, l) = fun_b(l/B^j, B);
    end
end
    
for j = 0:j_max
    disp(['j = ', num2str(j), ' starts...']);
    %compute Nside by the inequality Nside>=[B^{j+1}]/2
    Nside = 2^max((ceil(log2(fix(B^(j+1)))-1)), 0);
    disp(['Nside = ', num2str(Nside)])
    
    Nring = 4*Nside-1;
    pixList = inRing(Nside, 1:Nring);
    tp = pix2ang(Nside, 'nest', false);
    thetas = zeros(2*Nside,1);
    for r = 1:2*Nside
        thetas(r) = tp{pixList{r}(1)}(1);
    end

    pre_legendre = cell(l_max,1);
    for l = 1:l_max
        temp_mat = legendre(l, cos(thetas),'norm');
        temp_mat2 = fliplr(temp_mat(:,1:(2*Nside-1)));
        for m = 0:l
           temp_mat2(m+1,:) = (-1)^bitget(l+m, 1)*temp_mat2(m+1,:);
        end
        pre_legendre{l} = [temp_mat temp_mat2];
    end

    %compute the total number of cubature points Npix
    Npix = 12*Nside^2;
    %compute the cubature weights
    lambda = 4*pi/Npix;
    
    %compute the minimum index and the maximum index of l 
    %l_st>=1
    l_st = ceil(B^(j-1));
    l_en = min( floor(B^(j+1)), l_max);
    
    alm = coef;
    for l = l_st:l_en
        alm(l+1, bw:l+bw) = coef(l+1, bw:l+bw)*b_vector(j+1, l)*sqrt(lambda);
    end
    
    beta{j+1} = inv_spharmonic_tran( alm, Npix, Nring, pre_legendre, tp, pixList, bw, l_st, l_en );
   
end

e = cputime-t;
disp(['Elapsed CPU time: ', num2str(e)])

end
    
    
    