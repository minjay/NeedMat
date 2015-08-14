function A = get_A(B, j_min, j_max, theta, phi, n_dist)
%GET_A   Computes the design matrix A.
%
%   A = get_A(B, j_min, j_max, theta, phi, n_dist)
%
% Inputs:
%   B - the parameter
%   j_min - the minimal frequency
%   j_max - the maximal frequency
%   theta - the co-latitude of the locations, N-by-1 vector
%   phi - the longitude of the locations, N-by-1 vector
%   n_dist - the number of points on the fine grid
%
% Outputs:
%   A - the design matrix, N-by-M matrix
%
% Author: Minjie Fan, 2015

N = length(theta);
[x, y, z] = sph2cart(phi, pi/2-theta, 1);
l_max = floor(B^(j_max+1));
% column vector
dist = linspace(-1, 1, n_dist)';
P = p_polynomial_value( n_dist, l_max, dist );

bl_vector = get_bl_vector(B, j_max, l_max);

len_j = j_max-j_min+1;
psi = cell(len_j, 1);
Nside = zeros(len_j, 1);
Npix = zeros(len_j, 1);
for j = j_min:j_max
    Nside(j-j_min+1) = get_Nside(B, j);
    Npix(j-j_min+1) = 12*Nside(j-j_min+1)^2;
    sqrt_lambda = sqrt(4*pi/Npix(j-j_min+1));
    psi{j-j_min+1} = spneedlet_eval_fast(B, j, bl_vector, P, dist, sqrt_lambda);
end

M = sum(Npix);
A = zeros(N, M);
obs = zeros(N, 1);
index = 0;
for j = 1:len_j
    j
    tp = pix2ang(Nside(j), 'nest', false); 
    for k = 1:Npix(j)
        index = index+1;
        theta_xi = tp{k}(1);
        phi_xi = tp{k}(2);
        [x_xi, y_xi, z_xi] = sph2cart(phi_xi, pi/2-theta_xi, 1);
        for i = 1:N
            obs(i) = sum([x(i) y(i) z(i)].*[x_xi y_xi z_xi]);
            obs(i) = min(obs(i), 1);
            obs(i) = max(obs(i), -1);
        end
        A(:, index) = interp1(dist, psi{j}, obs, 'spline');
    end
end

end