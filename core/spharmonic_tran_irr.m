function alm = spharmonic_tran_irr(theta, phi, f, l_max)

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