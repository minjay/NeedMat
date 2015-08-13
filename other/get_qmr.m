function qmr = get_qmr(blmj, Nring, legendre_norm, bw, l_st, l_en)

qmr = zeros(l_en, Nring);
for m = 1:l_en
    l_st2 = max(m, l_st);
    legendre_mat = zeros(l_en-l_st2+1, Nring);
    for l = l_st2:l_en
        legendre_mat(l-l_st2+1, :) = legendre_norm{l}(m+1, :);
    end
    qmr(m, :) = reshape(blmj(l_st2+1:l_en+1, m+bw), 1, l_en-l_st2+1)*legendre_mat;
end

end