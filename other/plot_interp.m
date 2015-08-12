function plot_interp(theta_mat, phi_mat, f_mat, res_interp)

theta_interp = linspace(0, pi, res_interp/2);
phi_interp = linspace(0, 2*pi, res_interp);
[phi_mat_interp, theta_mat_interp] = meshgrid(phi_interp, theta_interp);

f_wend_hat_interp = interp2(phi_mat, theta_mat, f_mat, phi_mat_interp, theta_mat_interp, 'spline');
[L, T] = meshgrid(phi_interp-pi, pi/2-theta_interp);
[HX, HY] = sph2hammer(L, T);
figure
subplot('position',[0 0 1 1]);
pcolor(HX, HY, f_wend_hat_interp);
shading interp
axis equal
axis tight
axis off
colorbar('southoutside')

end