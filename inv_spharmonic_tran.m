function map = inv_spharmonic_tran( alm, Npix, Nring, pre_legendre, tp, pixList, bw, l_st, l_en )

%init
map = zeros(1, Npix);

temp_mat = zeros(l_en-l_st+1, Nring);
for l = l_st:l_en
    temp_mat(l-l_st+1,:) = pre_legendre{l}(1,:);
end
%No conj here!
term1 = alm(l_st+1:l_en+1, bw)'*temp_mat/sqrt(2*pi);

temp_mat2 = zeros(l_en, Nring);
for m = 1:l_en
    l_st2 = max(m, l_st);
    temp_mat = zeros(l_en-l_st2+1, Nring);
    for l = l_st2:l_en
        temp_mat(l-l_st2+1,:) = pre_legendre{l}(m+1,:);
    end
    %we need conj here!!!
    temp_mat2(m,:) = conj(alm(l_st2+1:l_en+1, m+bw)')*temp_mat/sqrt(2*pi)*(-1)^bitget(m, 1);
end

for r = 1:Nring
    for k = 1:length(pixList{r})
        index = pixList{r}(k);
        phi = tp{index}(2);
        temp_vec = exp((1:l_en)*1i*phi);
        product = temp_vec*temp_mat2(:,r);
        map(index) = real(term1(r))+2*real(product);
    end
end


    