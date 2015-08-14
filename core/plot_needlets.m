function plot_needlets(B, j, k, res)
%PLOT_NEEDLETS   Plots the spherical needlet with subscripts j and k.
%
%   plot_needlets(B, j, k, res)
%
% Inputs:
%   B - the parameter
%   j, k - the subscripts
%   res - the resolution of the plot
%
% Author: Minjie Fan, 2015

theta = linspace(0, pi, res/2);
phi = linspace(0, 2*pi, res);

[phi_mat, theta_mat] = meshgrid(phi, theta);
theta_vec = theta_mat(:);
phi_vec = phi_mat(:);

[dist, psi] = get_psi(B, j, k, theta_vec, phi_vec);

f = reshape(psi, res/2, res);

subplot('position',[0.1 0 0.85 1]);
[dist, I] =sort(acos(dist));
plot(dist, psi(I), 'LineWidth', 2)
axis equal
axis tight
xlabel('Great-circle distance', 'FontSize', 15)
ylabel('Amplitude', 'FontSize', 15)

[L, T] = meshgrid(phi-pi, pi/2-theta);
[HX, HY] = sph2hammer(L, T);
figure
subplot('position',[0 0 1 1]);
pcolor(HX, HY, f);
axis equal
axis tight
axis off
shading flat
colorbar('southoutside')

end
