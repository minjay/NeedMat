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

rng(1)
theta = theta+randn(n, 1)*pi/10;
theta(theta<0) = theta(theta<0)+pi;
theta(theta>pi) = theta(theta>pi)-pi;
phi = phi+randn(n, 1)*2*pi/10;
phi(phi<0) = phi(phi<0)+2*pi;
phi(phi>2*pi) = phi(phi>2*pi)-2*pi;

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
f_wend = (1-adj_dist).^6.*(35*adj_dist.^2+18*adj_dist+3)/3;
f_wend(adj_dist>1) = 0;

% plot
figure
subplot('position',[0 0 1 1]);
[HX, HY] = sph2hammer(phi-pi, pi/2-theta);
scatter(HX, HY, [], f_wend, 'filled')
colorbar('southoutside')
hold on

% draw contour
th = linspace(-pi/2,pi/2,101);
lam = -pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k', 'LineWidth', 2);
lam = pi+0*th;
[xh,yh] = sph2hammer(lam,th);
plot(xh,yh,'k', 'LineWidth', 2);

axis equal
axis off

l_max = 16;
alm = spharmonic_tran_irr(theta, phi, f_wend, l_max);

res = 100;

theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

[phi_mat, theta_mat] = meshgrid(phi, theta);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

f_wend_hat = inv_spharmonic_tran_naive(alm, theta_vec, phi_vec, l_max);
f_wend_hat = reshape(f_wend_hat, res/2, res);

res_interp = 1000;
plot_interp(theta_mat, phi_mat, f_wend_hat, res_interp)

% needlet transform
beta = spneedlet_tran(alm, l_max, B);

j_max = length(beta)-1;
n_j = zeros(j_max+1, 1);

for j = 0:j_max
    Nside = get_Nside(B, j);
    n_j(j+1) = 12*Nside^2;
end

n_dist = 1e3;
A = get_A(B, 0, j_max, theta_vec, phi_vec, n_dist);

f = cell(j_max+1, 1);
% hard thresholding
beta_trunc = beta;
q5 = quantile(abs(beta_trunc{5}), 0.95);
beta_trunc{5}(abs(beta_trunc{5})<q5) = 0;
for j = 0:j_max
    f{j+1} = A(:, sum(n_j(1:j+1))-n_j(j+1)+1:sum(n_j(1:j+1)))*beta_trunc{j+1};
end

f_wend_hat2 = alm(1, l_max+1)*sqrt(1/4/pi)+f{1}+f{2}+f{3}+f{4}+f{5};
f_wend_hat2 = reshape(f_wend_hat2, res/2, res);

plot_interp(theta_mat, phi_mat, f_wend_hat2, res_interp)
