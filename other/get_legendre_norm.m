function legendre_norm = get_legendre_norm(l_max, Nside) 

Nring = 4*Nside-1;
pixList = inRing(Nside, 1:Nring);
tp = pix2ang(Nside, 'nest', false);

theta = zeros(2*Nside, 1);
for r = 1:2*Nside
    theta(r) = tp{pixList{r}(1)}(1);
end

% \tilde P_lm
legendre_norm = cell(l_max, 1);

for l = 1:l_max
    Plm_tilde = legendre(l, cos(theta), 'norm');
    for m = 0:l
        Plm_tilde(m+1, :) = (-1)^bitget(m, 1)*Plm_tilde(m+1, :);
    end
    Plm_tilde2 = fliplr(Plm_tilde(:, 1:(2*Nside-1)));
    for m = 0:l
       Plm_tilde2(m+1,:) = (-1)^bitget(l+m, 1)*Plm_tilde2(m+1, :);
    end
    legendre_norm{l} = [Plm_tilde Plm_tilde2]/sqrt(2*pi);
end

end