function tautr = get_tautr(qmr, l_en, nr, phi0r)

tautr = zeros(nr, 1);

for m = 1:l_en
    t = mod(m, nr);
    tautr(t+1) = tautr(t+1)+qmr(m)*exp(1i*phi0r*m);
end

end
