clear;

% Channel with Orography Initial Conditions
% with approximately hydrostatic balance
% Vertical coefficients use Method 3 proposed by Li Yiyuan
% Mount profile from MAY WONG: Testing of a Cell-Integrated 
% Semi-Lagrangian Semi-Implicit Nonhydrostatic
% Atmospheric Solver (CSLAM-NH) with Idealized Orography

NX_full = 200; NX_half = NX_full + 1;
NLEV_full = 30; NLEV_half = NLEV_full + 1;

% constants
g = 9.80616; 
gamma = 1.4; Cp = 1004.0; Rd = 287.04;
% mount profile: h(x) = h0*[cos(pi*x/lambda)]^2*[cos(pi*x/2/a)]^2
h0 = 3.0e3/3*2; a = 25.0e3/3*2; lambda = 8.0e3/3*2;
% % mount profile: h(x) = h0*exp[-(x/a)^2]*cos[(pi*x/lambda)]^2
% h0 = 0.25e3; a = 5.0e3; lambda = 4.0e3;

% ref state
u_ref = 0.0; w_ref = 0.0; 
T_ref = 300.0; p_ref = 1.0e5;
N_freq = g / sqrt(Rd * T_ref);

% domain size:
x_lo = -150.0e3/3*2; x_hi = 150.0e3/3*2;
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
for i = 1:1:NX_half
    if abs(x_half(i)) < a
        zs_half(i) = h0 * (cos(pi*x_half(i)/lambda))^2 ...
                        * (cos(pi*x_half(i)/2/a))^2;
    else
        zs_half(i) = 0.0;
    end
end
for i = 1:1:NX_full
    if abs(x_full(i)) < a
        zs_full(i) = h0 * (cos(pi*x_full(i)/lambda))^2 ...
                        * (cos(pi*x_full(i)/2/a))^2;
    else
        zs_full(i) = 0.0;
    end
end
% zs_half = h0 * exp(-(x_half/a).^2) .* cos((pi*x_half/lambda)).^2;
% zs_full = h0 * exp(-(x_full/a).^2) .* cos((pi*x_full/lambda)).^2;

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
ph_top = ph_ns(1,1)

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
z_top = z(1,1)
gz = g * z;


% % 输出算例控制文件
% nstep = 180; nstep_to_write = 10;
% delta_t = 1.0;
% zh = 2./3. * z_top; tau_0 = 6.0;
% fid = fopen("../cmake-build-debug/case-control.txt", "w+");
% fprintf(fid, "method: 3\n");
% fprintf(fid, "x_lo: %f\n", x_lo);
% fprintf(fid, "x_hi: %f\n", x_hi);
% fprintf(fid, "NX_full: %d\n", NX_full);
% fprintf(fid, "NLEV_full: %d\n", NLEV_full);
% fprintf(fid, "nstep: %d\n", nstep);
% fprintf(fid, "nstep_to_write: %d\n", nstep_to_write);
% fprintf(fid, "delta_t: %f\n", delta_t);
% fprintf(fid, "u_ref: %f\n", u_ref);
% fprintf(fid, "w_ref: %f\n", w_ref);
% fprintf(fid, "Rayleigh_layer_height: %f\n", zh);
% fprintf(fid, "tau_0: %f\n", tau_0);
% fclose(fid);
% 
% % 输出初始场文件
% save("../cmake-build-debug/A_half.txt","A_half", "-ascii");
% save("../cmake-build-debug/B_half.txt","B_half", "-ascii");
% save("../cmake-build-debug/topo.txt","zs_half", "-ascii");
% save("../cmake-build-debug/phs0.txt","phs_full", "-ascii");
% save("../cmake-build-debug/u0.txt","u", "-ascii");
% save("../cmake-build-debug/w0.txt","w", "-ascii");
% save("../cmake-build-debug/pt0.txt","pt", "-ascii");
% save("../cmake-build-debug/geo0.txt","gz", "-ascii");

plot_width_ratio = 3;
subplot(1,2,1);
plot(x_full, gz); xlabel("x/m"); ylabel("gz/(m^2/s^2)");
xlim([-plot_width_ratio*a, plot_width_ratio*a]); title("geo potential"); 
subplot(1,2,2);
plot(x_full, ph_cell); xlabel("x/m"); ylabel("ph/Pa");
xlim([-plot_width_ratio*a, plot_width_ratio*a]); title("hydrostatic pressure"); 

