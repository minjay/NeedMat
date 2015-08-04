function Y_theta = spharmonic_partial_theta_eval( l, m, theta, phi )

signM = sign(m);
m = abs(m);

if m+1<=l
    tmp = spharmonic_eval(l, m+1, theta, phi);
else
    tmp = 0;
end
Y_theta = m*cot(theta)*spharmonic_eval(l, m, theta, phi)+...
    sqrt((l-m)*(l+m+1))*exp(-1i*phi)*tmp;

if 0>signM
    Y_theta = (-1)^m*conj(Y_theta);
end

end