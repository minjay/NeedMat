function beta = spneedlet_tran(coef, l_max, B)
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

%compute j_max from the inequality ceil(B^{j-1})<=l_max
j_max = get_j_max(B, l_max);

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
    Nside = 1;
    lb = floor(B^(j+1));
    while 2*Nside<lb
        Nside = Nside*2;
    end
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
    
    
    