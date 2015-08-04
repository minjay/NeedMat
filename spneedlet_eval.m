function psi = spneedlet_eval( theta, phi, B, j_max )
%SPNEEDLET_EVAL Needlet evaluation.
%
%   spneedlet_eval( theta, phi, B, j_max ) evaluates the spherical needlets psi_jk,
%   j=0,1,...,j_max, k=1,2,...,Npix, at the point (theta, phi) with the parameter B.
%   The return value psi is a cell.

psi = cell(j_max+1,1);

l_max = floor(B^(j_max+1));
bw = l_max+1;

%b_vector records the evaluations of the window function b
b_vector = zeros(j_max+1,l_max);
for j = 0:j_max
    for l = 1:l_max
        b_vector(j+1, l) = fun_b(l/B^j, B);
    end
end

coef = zeros(bw, 2*bw-1);

for l = 1:l_max
    for m = 0:l
        coef(l+1, m+bw) = conj(spharmonic_eval(l, m, theta, phi));
    end
end
    
for j = 0:j_max
    
    Nside = 2^max((ceil(log2(fix(B^(j+1)))-1)), 0);
    
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
    l_st = ceil(B^(j-1));
    l_en = floor(B^(j+1));
    
    alm = coef;
    
    for l = l_st:l_en
        alm(l+1, bw:l+bw) = coef(l+1, bw:l+bw)*b_vector(j+1, l)*sqrt(lambda);
    end
    
    psi{j+1} = inv_spharmonic_tran( alm, Npix, Nring, pre_legendre, tp, pixList, bw, l_st, l_en );
   
end

end