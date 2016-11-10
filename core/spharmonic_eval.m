function Y = spharmonic_eval(l, m, theta, phi)
%SPHARMONIC_EVAL   Evaluates the spherical harmonic with subscripts l and m 
%at locations (theta, phi).
%
%   Y = spharmonic_eval(l, m, theta, phi)
%
% Inputs:
%   l, m - the subscripts
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%
% Outputs:
%   Y - the values of the spherical harmonic at the locations, N-by-1 vector
%
% Author: Minjie Fan, 2015

signM = sign(m);
m = abs(m);
r = bitget(m, 1);

C = (-1)^r/sqrt(2*pi);

P = legendre( l, cos(theta), 'norm' );

Y = C*P(m+1, :)'.*exp(1i*m*phi);

if 0>signM
    Y = (-1)^r*conj(Y);
end

end