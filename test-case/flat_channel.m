clear;

% Flat Channel Initial Conditions
% Steady State with approximately hydrostatic balance 
% with no Orography

NX_full = 10; NX_half = NX_full + 1;
NLEV_full = 12; NLEV_half = NLEV_full + 1;
% eta = ((ph-ph_top)/(phs-ph_top))^(1/n)
n_power = 1.1;
% constants
g = 9.80616; 
gamma = 1.4; Cp = 1004.0; Rd = 287.04;
% mount profile: flat channel
% ref state
u_ref = 0.0; w_ref = 0.0; 
T_ref = 300.0; p_ref = 1.0e5;
N_freq = g / sqrt(Rd * T_ref);
% domain size
x_lo = -5.0e3; x_hi = 5.0e3;
z_top = 12.0e3;
ph_top = 205.448e2;
diss_rate = 0.006511; % 温度衰减率 K/m

% 计算水平坐标
delta_X = (x_hi-x_lo)/NX_full;
x_half = x_lo : delta_X : x_hi;
for i = 1:1:NX_full
    x_full(i) = (x_half(i)+x_half(i+1))/ 2.0;
end

% 计算垂直坐标
delta_eta = 1.0 / NLEV_full;
eta_half = 0.0 : delta_eta : 1.0;
eta_half = transpose(eta_half);
for k = 1:NLEV_full
    eta_full(k) = 0.5 * ( eta_half(k) + eta_half(k+1) );
end
eta_full = transpose(eta_full);
B_half = eta_half.^n_power;
B_full = eta_full.^n_power;
A_half = (1.0 - B_half) / p_ref * ph_top; 
A_full = (1.0 - B_full) / p_ref * ph_top;

% 计算地表高度
zs_half = zeros(1,NX_half);
zs_full = zeros(1,NX_full);

% 计算地面静力气压
phs_full = p_ref * (1.0 - diss_rate/T_ref * zs_full).^(g/Rd/diss_rate);

% 计算各个位置full level的静力气压
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        ph_cell(k,i) = A_full(k)*p_ref + B_full(k)*phs_full(i);
    end
end
% 计算各个位置half level的静力气压
for k = 1:1:NLEV_half
    for i = 1:1:NX_full
        ph_ns(k,i) = A_half(k)*p_ref + B_half(k)*phs_full(i);
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
gz = g * z;

% 输出算例控制文件
nstep = 48000; nstep_to_write = 2400;
delta_t = 1.5;
zh = 2./3. * z_top; tau_0 = 6.0;
fid = fopen("../cmake-build-debug/case-control.txt", "w+");
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
save("../cmake-build-debug/A_full.txt","A_full", "-ascii");
save("../cmake-build-debug/B_full.txt","B_full", "-ascii");
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

