clear

% 画图用的

w = load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/w_36000.00s.dat");
u = load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/u_36000.00s.dat") - 10;
gz= load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/gz_36000.00s.dat");
u_ns = load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/u_ns_36000.00s.dat") - 10;

x_full = load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/x_full.dat");
eta_full =load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/eta_full.dat");
x_half= load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/x_half.dat");
eta_half =load("../cmake-build-release/schar_10h_1.0s_K不迎风_位温1.0_水平500m_隐式R耗散_无水平和三维散度耗散/eta_half.dat");

NLEV_full = size(eta_full); NLEV_full = NLEV_full(1); 
NLEV_half = size(eta_half); NLEV_half = NLEV_half(1); 
NX_full = size(x_full); NX_full = NX_full(2); 
NX_half = size(x_half); NX_half = NX_half(2);

delta_z = 50; % m
z_phys = delta_z:delta_z:12.0e3;
num_interp_level = size(z_phys); num_interp_level = num_interp_level(2);

x_ori = repmat(x_full, NLEV_half, 1);
z_ori = gz / 9.80616;
delta_x = 250;
x_phys = -200e3+delta_x/2:delta_x:200e3-delta_x/2;
num_interp_hori = size(x_phys); num_interp_hori = num_interp_hori(2);
x_interp = repmat(x_phys, num_interp_level, 1);
z_interp = repmat(z_phys', 1, num_interp_hori);

w_interp = griddata(x_ori,z_ori,w, x_interp,z_interp);
u_ns_interp = griddata(x_ori,z_ori,u_ns, x_interp,z_interp);

subplot(1,2,1);
w_value_start = -1;
w_value_end   = 1;
w_values = w_value_start:0.05:w_value_end;
contour(x_interp/1e3, z_interp/1e3, w_interp, w_values, 'LineWidth',1.5);
caxis([w_value_start, w_value_end]);
set(gca,'xlim',[-20,20]); set(gca,'ylim',[0,12]);
set(gca,'XTick',[-20,-10,0,10,20]);
xlabel("Horizontal distance (km)");
ylabel("Height (km)");
title("Vertical velocity: Schar mountain wave after T = 10h");

subplot(1,2,2);
u_value_start = -2;
u_value_end   = 2;
u_values = u_value_start:0.2:u_value_end;
contour(x_interp/1e3, z_interp/1e3, u_ns_interp, u_values, 'LineWidth',1.5);
caxis([u_value_start, u_value_end]);
set(gca,'xlim',[-20,20]); set(gca,'ylim',[0,12]);
set(gca,'XTick',[-20,-10,0,10,20]);
xlabel("Horizontal distance (km)");
ylabel("Height (km)");
title("Horizontal velocity: Schar mountain wave after T = 10h");
