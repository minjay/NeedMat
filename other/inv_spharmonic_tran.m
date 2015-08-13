function map = inv_spharmonic_tran(blmj, Nside, legendre_norm, bw, l_st, l_en)

% init
Npix = 12*Nside^2;
Nring = 4*Nside-1;
pixList = inRing(Nside, 1:Nring);
tp = pix2ang(Nside, 'nest', false);
n = zeros(Nring, 1);
phi0 = zeros(Nring, 1);
for r = 1:Nring
    n(r) = length(pixList{r});
    phi0(r) = tp{pixList{r}(1)}(2);
end
map = zeros(Npix, 1);

legendre_mat = zeros(l_en-l_st+1, Nring);
for l = l_st:l_en
    legendre_mat(l-l_st+1,:) = legendre_norm{l}(1,:);
end
term1 = reshape(blmj(l_st+1:l_en+1, bw), 1, l_en-l_st+1)*legendre_mat;

qmr = get_qmr(blmj, Nring, legendre_norm, bw, l_st, l_en);

index = 0;
for r = 1:Nring
    tautr = get_tautr(qmr(:, r), l_en, n(r), phi0(r));
    term2 = ifft(tautr)*n(r); 
    for p = 1:n(r)
        index = index+1;
        map(index) = term1(r)+2*real(term2(p));
    end
end

end
    