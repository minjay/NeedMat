function Y = spharmonic_eval( l, m, theta, phi )
%SPHARMONIC_EVAL Evaluates the value of the spherical harmonic Y_{lm} at
%the point (theta, phi)
%
%   spharmonic_eval( l, m, theta, phi ) evaluates the value of the
%spherical harmonic Y_{lm} at the point (theta, phi), where l is the
%degree and m is the order. Note that Y_{lm} is expressed in the complex
%form, and there is Condon-Shortley phase (-1)^m in the associated Legendre
%polynomials. Besides, the angular unit measuring theta and phi is radian.
%
%   This function depends on the Matlab built-in function legendre.
%
%   Written by Minjie Fan

signM = sign(m);
m = abs(m);

C = (-1)^bitget(m, 1)/sqrt(2*pi);

P = legendre( l, cos(theta), 'norm' );

Y = C*P(m+1)*exp(1i*m*phi);

if 0>signM
    Y = (-1)^m*conj(Y);
end

end