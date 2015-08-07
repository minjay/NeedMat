function plot_needlets(B, j, res)

theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

j_max = j;
l_max = floor(B^(j_max+1));

bl_vector = get_bl_vector(B, j_max, l_max);

Nside = 1;
lb = floor(B^(j+1));
while 2*Nside<lb
    Nside = Nside*2;
end

sqrt_lambda = sqrt(4*pi/(12*Nside^2));
tp = pix2ang(Nside, 'nest', false);

k = length(tp)/2;

theta_xi = tp{k}(1);
phi_xi = tp{k}(2);
[x_xi, y_xi, z_xi] = sph2cart(phi_xi, pi/2-theta_xi, 1);

[phi_vec, theta_vec] = meshgrid(phi, theta);
theta_vec = theta_vec(:);
phi_vec = phi_vec(:);
[x, y, z] = sph2cart(phi_vec, pi/2-theta_vec, 1);

n = length(x);
dist = zeros(n, 1);

for i = 1:n
    dist(i) = sum([x(i) y(i) z(i)].*[x_xi y_xi z_xi]);
end

P = p_polynomial_value( n, l_max, dist );

psi = spneedlet_eval_fast(B, j, bl_vector, P, dist, sqrt_lambda);

f = reshape(psi, res/2, res);

figure
plot(acos(dist(:)), f(:), '.')

theta = pi/2-theta;
phi = phi-pi;
[L, T] = meshgrid(phi,theta);
[HX, HY] = sph2hammer(L, T);
figure
pcolor(HX, HY, f);
axis equal
axis tight
shading flat
colorbar

end
