clear;

% Channel with Orography Initial Conditions
% with approximately hydrostatic balance
% Vertical coefficients use Method 3 proposed by Li Yiyuan
% from Section 3.a, May Wong: Testing of a Cell-Integrated
% Semi-Lagrangian Semi-Implicit Nonhydrostatic Atmospheric 
% Solver (CSLAM-NH) with Idealized Orography

NX_full = 400; NX_half = NX_full + 1;
NLEV_full = 40; NLEV_half = NLEV_full + 1;

% constants
g = 9.80616; 
gamma = 1.4; Cp = 1004.0; Rd = 287.04;
% mount profile: h = hc * ac^2 / (ac^2 + (x-xc)^2)
hc = 1000.0; xc = 0.0; ac = 2.0e3;

% ref state
u_ref = 0.0; w_ref = 0.0; 
T_ref = 300.0; p_ref = 1.0e5;
N_freq = g / sqrt(Rd * T_ref);

% domain size:
x_lo = -72.0e3; x_hi = 72.0e3;
ph_top_set = 0.0;
diss_rate = 0.006511; % 温度衰减率 K/m

% 计算水平坐标
delta_X = (x_hi-x_lo)/NX_full;
x_half = x_lo : delta_X : x_hi;
for i = 1:1:NX_full
    x_full(i) = (x_half(i)+x_half(i+1))/ 2.0;
end

% 计算垂直坐标 
eta_c = 0.2; % 模式顶部的eta：对于30km模式顶一般取0.2
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
% num_full = 1:1:NLEV_full; num_half = 1:1:NLEV_half;
% plot(num_half, eta_half, num_half, A_half, num_half, B_half);
% legend("eta_{half}","A_{half}","B_{half}");
% plot(num_full, eta_full, num_full, A_full, num_full, B_full);
% legend("eta_{full}","A_{full}","B_{full}");

% 计算地表高度
zs_half = hc * ac^2 ./ (ac^2 + (x_half-xc).^2);
zs_full = hc * ac^2 ./ (ac^2 + (x_full-xc).^2);
% plot(x_half, zs_half);

% 计算地面静力气压
phs_full = p_ref * (1.0 - diss_rate/T_ref * zs_full).^(g/Rd/diss_rate);

% 计算各个位置full level的静力气压
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        ph_cell(k,i) = A_full(k)*(p_ref-ph_top_set) ...
                + B_full(k)*(phs_full(i)-ph_top_set) + ph_top_set;
    end
end
% 计算各个位置half level的静力气压
for k = 1:1:NLEV_half
    for i = 1:1:NX_full
        ph_ns(k,i) = A_half(k)*(p_ref-ph_top_set) ...
                + B_half(k)*(phs_full(i)-ph_top_set) + ph_top_set;
    end
end

% 计算水平速度
u = zeros(NLEV_full, NX_half);
for k = 1:1:NLEV_full
    for i = 1:1:NX_half
        u(k,i) = u_ref;
    end
end

% 计算垂直速度
w = zeros(NLEV_half, NX_full);
for k = 1:1:NLEV_half
    for i = 1:1:NX_full
        w(k,i) = w_ref;
    end
end

% 计算位温
T = T_ref * (ph_cell / p_ref).^(Rd*diss_rate/g);
pt = T .* power(ph_cell / p_ref, (1-gamma) / gamma);

% 计算地势
z = T_ref/diss_rate * ( 1.0 - (ph_ns/p_ref).^(Rd*diss_rate/g) );
% 注意在很高的高空中，温度衰减率可能已经不是这个数了
gz = g * z;

% determine where to truncate
z_top_set = 12.0e3; 
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

% 输出算例控制文件
nstep = 18000; nstep_to_write = 300;
delta_t = 0.2;
zh = 1./2. * z_top_calc; tau_0 = 6.0;
fid = fopen("../cmake-build-debug/case-control.txt", "w+");
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
fclose(fid);

% 输出初始场文件
save("../cmake-build-debug/A_half.txt","A_half", "-ascii");
save("../cmake-build-debug/B_half.txt","B_half", "-ascii");
save("../cmake-build-debug/topo.txt","zs_half", "-ascii");
save("../cmake-build-debug/phs0.txt","phs_full", "-ascii");
save("../cmake-build-debug/u0.txt","u", "-ascii");
save("../cmake-build-debug/w0.txt","w", "-ascii");
save("../cmake-build-debug/pt0.txt","pt", "-ascii");
save("../cmake-build-debug/geo0.txt","gz", "-ascii");

subplot(1,2,1);
plot(x_full, gz); xlabel("x/m"); ylabel("gz/(m^2/s^2)");
title("geo potential"); 
subplot(1,2,2);
plot(x_full, ph_cell); xlabel("x/m"); ylabel("ph/Pa");
title("hydrostatic pressure"); 

