function alm = spharmonic_tran_irr(theta, phi, f, l_max)
%SPHARMONIC_TRAN_IRR   Spherical harmonic transform for irregularly spaced 
%observations of function f at locations (theta, phi). It estimates the 
%spherical harmonic coefficients empirically using the weights determined 
%by the Voronoi diagram of these locations on the sphere.
%
%   alm = spharmonic_tran_irr(theta, phi, f, l_max)
%
% Inputs:
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%   f - the values of the function at the locations
%   l_max - the maximal l
%
% Outputs:
%   alm - the spherical harmonic coefficients, (l_max+1)-by-(2*l_max+1)
%   matrix, alm(l+1, m+l_max+1) is the spherical harmonic coefficient with
%   subscripts l and m
%
% Author: Minjie Fan, 2015

n = length(theta);

% get weights
dirs = [phi, pi/2-theta];
faces = sphDelaunay(dirs);
[voronoi, duplicates] = sphVoronoi(dirs, faces);
voronoi = sphVoronoiAreas(voronoi);
area = voronoi.area;

% estimate alm
alm = zeros(l_max+1, 2*l_max+1);
tmp = zeros(n, 1);
for l = 0:l_max
    for m = 0:l
        for i = 1:n
            tmp(i) = conj(spharmonic_eval(l, m, theta(i), phi(i)));
        end
        alm(l+1, m+l_max+1) = sum(area.*f.*tmp);
    end
end

end