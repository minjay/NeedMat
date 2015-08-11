function Y = spharmonic_eval( l, m, theta, phi )
%SPHARMONIC_EVAL Evaluates the value of the spherical harmonic Y_{lm} at
%the point (theta, phi)

signM = sign(m);
m = abs(m);
r = bitget(m, 1);

C = (-1)^r/sqrt(2*pi);

P = legendre( l, cos(theta), 'norm' );

Y = C*P(m+1)*exp(1i*m*phi);

if 0>signM
    Y = (-1)^r*conj(Y);
end

end