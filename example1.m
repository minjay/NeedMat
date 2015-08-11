clear

% the grid
B = 2;
Nside = 8;
tp = pix2ang(Nside, 'nest', false);
n = length(tp);
theta = zeros(n, 1);
phi = zeros(n, 1);
for i = 1:n
    theta(i) = tp{i}(1);
    phi(i) = tp{i}(2);
end

theta = theta+randn(n, 1)*pi/10;
phi = phi+randn(n, 1)*2*pi/10;

% the center
theta_xi = pi/2;
phi_xi = pi;
[x_xi, y_xi, z_xi] = sph2cart(phi_xi, pi/2-theta_xi, 1);

[x, y, z] = sph2cart(phi, pi/2-theta, 1);

% the distance vector
dist = zeros(n, 1);

for i = 1:n
    dist(i) = acos(sum([x(i) y(i) z(i)].*[x_xi y_xi z_xi]));
end

% the support is pi/4
adj_dist = dist/(pi/4);
% the Wendland radial basis function
psi = (1-adj_dist).^6.*(35*adj_dist.^2+18*adj_dist+3)/3;
psi(adj_dist>1) = 0;

% plot
worldmap([-90 90], [-180 180])
scatterm((pi/2-theta)/pi*180, (phi-pi)/pi*180, [], psi, 'filled')

l_max = 16;
alm = spharmonic_tran_irr(theta, phi, l_max, psi);

psi_hat = inv_spharmonic_tran_irr(alm, theta, phi, l_max);

% needlet transform
beta = spneedlet_tran(alm, l_max, B);

j_max = length(beta)-1;
n_j = zeros(j_max+1, 1);

for j = 0:j_max
    Nside = get_Nside(B, j);
    tp = pix2ang(Nside, 'nest', false);
    n_j(j+1) = length(tp);
    theta_j = zeros(n_j(j+1), 1);
    phi_j = zeros(n_j(j+1), 1);
    for i = 1:n_j(j+1)
        theta_j(i) = tp{i}(1);
        phi_j(i) = tp{i}(2);
    end
    subplot(3, 2, j+1)
    worldmap([-90 90], [-180 180])
    scatterm((pi/2-theta_j)/pi*180, (phi_j-pi)/pi*180, 50-10*j, beta{j+1}, 'filled')
end

n_dist = 1e3;
A = get_A(B, 0, j_max, theta, phi, n_dist);

f = cell(j_max+1, 1);
beta_trunc = beta;
q5 = quantile(abs(beta_trunc{5}), 0.95);
beta_trunc{5}(abs(beta_trunc{5})<q5) = 0;
for j = 0:j_max
    f{j+1} = A(:, sum(n_j(1:j+1))-n_j(j+1)+1:sum(n_j(1:j+1)))*beta_trunc{j+1};
end

psi_hat2 = alm(1, l_max+1)*sqrt(1/4/pi)+f{1}+f{2}+f{3}+f{4}+f{5};

plot(dist, psi, '.')
hold on 
plot(dist, psi_hat, 'r.')
plot(dist, psi_hat2, 'g.')
