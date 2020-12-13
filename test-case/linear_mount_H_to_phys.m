clear

w = load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/w_64800.00s.dat");
u = load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/u_64800.00s.dat") - 10;
gz= load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/gz_64800.00s.dat");

x_full = load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/x_full.dat");
eta_full =load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/eta_full.dat");
x_half= load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/x_half.dat");
eta_half =load("../cmake-build-release/linear_H_18h_2.5s_位温1.0_隐式R耗散_无水平和三维散度耗散/eta_half.dat");

NLEV_full = size(eta_full); NLEV_full = NLEV_full(1); 
NLEV_half = size(eta_half); NLEV_half = NLEV_half(1); 
NX_full = size(x_full); NX_full = NX_full(2); 
NX_half = size(x_half); NX_half = NX_half(2);

delta_z = 50; % m
z_phys = delta_z:delta_z:12.0e3;
num_interp_level = size(z_phys); num_interp_level = num_interp_level(2);

x_ori = repmat(x_full, NLEV_half, 1);
z_ori = gz / 9.80616;
x_interp = repmat(x_full, num_interp_level, 1);
z_interp = repmat(z_phys', 1, NX_full);

w_interp = griddata(x_ori,z_ori,w, x_interp,z_interp);
% u_ns_interp = griddata(repmat(x,NLEV_half,1),z_ori,u_ns, x_interp,z_interp);
% values = [-0.015, -0.013, -0.011, -0.009, -0.007, -0.005, -0.003, -0.001, 0,...
%            0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015];
values = -0.016:0.002:0.016;
contourf(x_interp/1e3, z_interp/1e3, w_interp, values);
caxis([-0.015, 0.015]);
%set(gca,'XTick',-40e3:10e3:40e3); 
set(gca,'xlim',[-40, 40]); set(gca,'ylim',[0, 12]);
xlabel("Horizontal distance (km)");
ylabel("Height (km)");
title("Vertical velocity: Linear hydrostatic wave after T = 18h");

% subplot(1,2,1);
% contour(u_ns_interp,-2:0.2:2); title("u:250m pt 1.0"); 
% xlim([800-40,800+40]);
% subplot(1,2,2); 
% contour(w_interp,-2:0.05:2); title("w:250m pt 1.0");
% xlim([800-40,800+40]);