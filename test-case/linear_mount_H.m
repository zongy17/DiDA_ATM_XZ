clear;

% Linear Mountain Initial Conditions 
% from Section 3.a, May Wong: Testing of a Cell-Integrated
% Semi-Lagrangian Semi-Implicit Nonhydrostatic Atmospheric 
% Solver (CSLAM-NH) with Idealized Orography

NX_full = 120; NLEV_full = 100; % hydrostatic
NX_half = NX_full + 1; NLEV_half = NLEV_full + 1;

% constants
g = 9.80616; 
gamma = 1.4; Cp = 1004.0; Rd = 287.04;

% mount profile: linear mountain, unit: m
% h = hc * ac^2 / (ac^2 + (x-xc)^2)
hc = 10.0; xc = 0.0; ac = 10.0e3; % hydrostatic

% ref state
u_ref = 10.0; w_ref = 0.0; 
T_ref = 300.0; p_ref = 1.0e5;
N_freq = g / sqrt(Rd * T_ref);
% domain size
x_lo = -60.0e3; x_hi = 60.0e3; % hydrostatic
ph_top_set = 0.0;
diss_rate = 0.006511; % �¶�˥���� K/m

% ����ˮƽ����
delta_X = (x_hi-x_lo)/NX_full;
x_half = x_lo : delta_X : x_hi;
for i = 1:1:NX_full
    x_full(i) = (x_half(i)+x_half(i+1))/ 2.0;
end

% ���㴹ֱ���� 
eta_c = 0.1; % ģʽ������eta������30kmģʽ��һ��ȡ0.2
c1 = 2.0 * eta_c^2 / (1.0 - eta_c)^3;
c2 = - eta_c * (4.0 + eta_c + eta_c^2) / (1.0 - eta_c)^3;
c3 = 2.0 * (1.0 + eta_c + eta_c^2) / (1.0 - eta_c)^3;
c4 = - (1.0 + eta_c) / (1.0 - eta_c)^3;
eta_half = linspace(0.0, 1.0, NLEV_half);
eta_half = transpose(eta_half);
for k = 1:1:NLEV_half
    if eta_half(k,1) < eta_c
        B_half(k,1) = 0.0;
    else
        B_half(k,1) = c1 + c2 * eta_half(k,1)...
                + c3 * eta_half(k,1)^2 + c4 * eta_half(k,1)^3;
    end
end
A_half = eta_half - B_half;
for k = 1:1:NLEV_full
    A_full(k,1) = 0.5*( A_half(k,1)+A_half(k+1,1) );
    B_full(k,1) = 0.5*( B_half(k,1)+B_half(k+1,1) );
end
eta_full = A_full + B_full;

% ����ر�߶�
zs_half = hc ./ (1 + (x_half-xc)/ac .* (x_half-xc)/ac);
zs_full = hc ./ (1 + (x_full-xc)/ac .* (x_full-xc)/ac);
right = zs_half(1,2:NX_half);
left  = zs_half(1,1:NX_half-1);
dzsdx_full = (right-left) / delta_X;

% ������澲����ѹ
phs_full = p_ref * exp(- g * zs_full / Rd / T_ref);
% �������λ��full level�ľ�����ѹ
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        ph_cell(k,i) = A_full(k)*p_ref + B_full(k)*phs_full(i);
    end
end
% �������λ��half level�ľ�����ѹ
for k = 1:1:NLEV_half
    for i = 1:1:NX_full
        ph_ns(k,i) = A_half(k)*p_ref + B_half(k)*phs_full(i);
    end
end

% ����ˮƽ�ٶ�
u = zeros(NLEV_full, NX_half);
for k = 1:1:NLEV_full
    for i = 1:1:NX_half
        u(k,i) = u_ref;
    end
end

% ���㴹ֱ�ٶ�
w = zeros(NLEV_half, NX_full);
for k = 1:1:NLEV_half
    for i = 1:1:NX_full
        w(k,i) = w_ref;
    end
end

% ����λ��
T = T_ref * (ph_cell / p_ref).^(Rd*diss_rate/g);
pt = T .* power(ph_cell / p_ref, (1-gamma) / gamma);

% �������
z = T_ref/diss_rate * ( 1.0 - (ph_ns/p_ref).^(Rd*diss_rate/g) );
gz = g * z;

% determine where to truncate
z_top_set = 20.0e3; % hydrostatic

trunc_lev = NLEV_half;
for k = 1:1:NLEV_half
    if z(k,1) < z_top_set
        trunc_lev = k;
        break
    end
end
% do truncation
A_half = A_half(trunc_lev:NLEV_half, :); 
B_half = B_half(trunc_lev:NLEV_half, :);
eta_half = eta_half(trunc_lev:NLEV_half, :);
A_full = A_full(trunc_lev:NLEV_full, :);
B_full = B_full(trunc_lev:NLEV_full, :);
eta_full = eta_full(trunc_lev:NLEV_full, :);
ph_ns = ph_ns(trunc_lev:NLEV_half, :);
ph_cell = ph_cell(trunc_lev:NLEV_full, :);
u = u(trunc_lev:NLEV_full, :);
w = w(trunc_lev:NLEV_half, :);
T = T(trunc_lev:NLEV_full, :);
pt = pt(trunc_lev:NLEV_full, :);
z  = z(trunc_lev:NLEV_half, :);
gz = gz(trunc_lev:NLEV_half, :);
NLEV_half = size(eta_half,1); NLEV_full = NLEV_half - 1;
z_top_calc = z(1,1) %#ok<NOPTS>
ph_top_calc = ph_ns(1,1) %#ok<NOPTS>

% determine the maximum delta_x and delta_z: for the time step
min_delta_x = inf; min_delta_z = inf;
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        local_delta_x = x_half(i+1) - x_half(i);
        local_delta_z = z(k,i) - z(k+1,i);
        if (local_delta_x < min_delta_x)
            min_delta_x = local_delta_x;
        end
        if (local_delta_z < min_delta_z)
            min_delta_z = local_delta_z;
        end
    end
end
min_delta_x %#ok<NOPTS>
min_delta_z %#ok<NOPTS>

% ������������ļ�
simu_time = 18.0;% Сʱ
delta_t = 2.5; % ��
delta_t_to_write = 30; % ����
nstep = simu_time * 3600 / delta_t; 
nstep_to_write = delta_t_to_write * 60 / delta_t;

kexi = 0.6; % kexi=0.5 for Crank-Nicholson
zh = z_top_calc - 10.0e3; tau_0 = 10.0; % hydrostatic

fid = fopen("../cmake-build-release/case-control.txt", "w+");
fprintf(fid, "method: 3\n");
fprintf(fid, "x_lo: %f\n", x_lo);
fprintf(fid, "x_hi: %f\n", x_hi);
fprintf(fid, "NX_full: %d\n", NX_full);
fprintf(fid, "NLEV_full: %d\n", NLEV_full);
fprintf(fid, "nstep: %d\n", nstep);
fprintf(fid, "nstep_to_write: %d\n", nstep_to_write);
fprintf(fid, "delta_t: %f\n", delta_t);
fprintf(fid, "u_ref: %f\n", u_ref);
fprintf(fid, "w_ref: %f\n", w_ref);
fprintf(fid, "Rayleigh_layer_height: %f\n", zh);
fprintf(fid, "tau_0: %f\n", tau_0);
fprintf(fid, "kexi: %f\n", kexi);
fclose(fid);

% �����ʼ���ļ�
save("../cmake-build-release/A_half.txt","A_half", "-ascii");
save("../cmake-build-release/B_half.txt","B_half", "-ascii");
save("../cmake-build-release/zs_half.txt","zs_half", "-ascii");
save("../cmake-build-release/zs_full.txt","zs_full", "-ascii");
save("../cmake-build-release/phs0.txt","phs_full", "-ascii");
save("../cmake-build-release/u0.txt","u", "-ascii");
save("../cmake-build-release/w0.txt","w", "-ascii");
save("../cmake-build-release/pt0.txt","pt", "-ascii");
save("../cmake-build-release/geo0.txt","gz", "-ascii");

subplot(1,2,1);
plot(x_full, gz); xlabel("x/m"); ylabel("gz/(m^2/s^2)");
title("geo potential"); 
subplot(1,2,2);
plot(x_full, ph_cell); xlabel("x/m"); ylabel("ph/Pa");
title("hydrostatic pressure"); 