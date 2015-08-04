function Y_phi = spharmonic_partial_phi_eval( l, m, theta, phi )

signM = sign(m);
m = abs(m);

Y_phi = 1i*m*spharmonic_eval(l, m, theta, phi);

if 0>signM
    Y_phi = (-1)^m*conj(Y_phi);
end

end