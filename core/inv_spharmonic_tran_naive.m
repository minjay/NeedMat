function map = inv_spharmonic_tran_naive(alm, theta, phi, l_max)

n = length(theta);
map = zeros(n, 1);
for i = 1:n
    if mod(i, 1000)==0
        i
    end
    for l = 0:l_max
        for m = 0:l
            if m==0
                map(i) = map(i)+alm(l+1, m+l_max+1)*spharmonic_eval(l, m, theta(i), phi(i));
            else
                map(i) = map(i)+2*real(alm(l+1, m+l_max+1)*spharmonic_eval(l, m, theta(i), phi(i)));
            end
        end
    end
end

end
