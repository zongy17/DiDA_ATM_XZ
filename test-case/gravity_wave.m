clear;

% Gravity Wave Initial Conditions 
% from Melvin, T., et al., An inherently mass-conserving 
% iterative semi-implicit semi-Lagrangian discretization of 
% the non-hydrostatic vertical-slice equations. 
% Quarterly journal of the Royal Meteorological Society, 2010. 136(648): p. 799-n/a.

% constants
g = 9.80616; 
gamma = 1.4; Cp = 1004.0; Rd = 287.04;
kappa = Rd/Cp;

% ref state
u_ref = 20.0; w_ref = 0.0; % m/s
p_ref = 1.0e5; % Pa 
pt_ref = 300; % K
T_surf = 300; % K, on the surface, pt = T
N_freq = 0.01; % s^(-1)

% domain size
% x_lo = -240.0e3; x_hi = 240.0e3; NX_full = 1920;
x_lo = -150e3; x_hi = 150e3; NX_full = 300;
ph_top_set = 0.0;

% 计算水平坐标
NX_half = NX_full + 1;
delta_X = (x_hi-x_lo)/NX_full;
x_half = x_lo : delta_X : x_hi;
for i = 1:1:NX_full
    x_full(i) = (x_half(i)+x_half(i+1))/ 2.0;
end

% 垂直坐标系数
NLEV_full = 10; NLEV_half = NLEV_full + 1;
if NLEV_full == 10
    A_half = [
          0.2737592327457158 % 1
          0.2574691274856134 % 2
          0.2392257778961169 % 3
          0.2188579122253837 % 4
          0.1961834482057335 % 5
          0.171008960410953  % 6
          0.1431291236866448 % 7
          0.1123261317013542 % 8
          0.07836908962461686 % 9
          0.04101337989532572 % 10
          0.0                 % 11
    ];
    B_half = [
          0.0                 % 1
          0.05950522689853417 % 2
          0.1261453522616924  % 3
          0.2005460052239693  % 4
          0.283372304056826   % 5
          0.3753308018298088  % 6
          0.4771715194731282  % 7
          0.5896900697201707  % 8
          0.7137298755603581  % 9
          0.8501844869888955  % 10
          1.0                 % 11
    ];
elseif NLEV_full == 40
    A_half = [
          0.1025039799580591 %  1
          0.1085106112775176 %  2
          0.1148692252207011 %  3
          0.1216004476194305 %  4
          0.1287261129587656 %  5
          0.1362693352028837 %  6
          0.1442545827712838 %  7
          0.1527077579085218 %  8
          0.1616562807049334 %  9
          0.1711291780408878 % 10
          0.1811571777430863 % 11
          0.1917728082583264 % 12
          0.1904697560360646 % 13
          0.1992215664069182 % 14
          0.2081033503253796 % 15
          0.2170765373305924 % 16
          0.2260947251739    % 17
          0.2351025487717675 % 18
          0.2440344018023399 % 19
          0.2528129925300951 % 20
          0.2613477131840651 % 21
          0.2695327996814114 % 22
          0.2772452556464864 % 23
          0.2843425114885985 % 24
          0.2906597857277424 % 25
          0.2960071117496997 % 26
          0.3001659886776721 % 27
          0.302885610008119  % 28
          0.3038786180077924 % 29
          0.3028163255331611 % 30
          0.299323339829651  % 31
          0.2929715149035437 % 32
          0.2832731501298948 % 33
          0.2696733427488229 % 34
          0.2515413906792607 % 35
          0.2281611294972407 % 36
          0.1987200733208144 % 37
          0.1622972135317635 % 38
          0.1178493115391018 % 39
          0.06419550191988976 % 40
          0.0                 % 41
    ];
    B_half = [
          0.0 %  1
          0.0 %  2
          0.0 %  3
          0.0 %  4
          0.0 %  5
          0.0 %  6
          0.0 %  7
          0.0 %  8
          0.0 %  9
          0.0 % 10
          0.0 % 11
          0.0 % 12
          0.01254074813198438 % 13
          0.01568515147901569 % 14
          0.01939668757548295 % 15
          0.02375477661820929 % 16
          0.02884906434573413 % 17
          0.03478069335699107 % 18
          0.0416637300061256  % 19
          0.04962676577066524 % 20
          0.05881471427746436 % 21
          0.06939082773441382 % 22
          0.08153895938918876 % 23
          0.09546610185670711 % 24
          0.111405234766373   % 25
          0.1296185192255581 % 26
          0.1504008811297395 % 27
          0.1740840304312517 % 28
          0.2010409691738044 % 29
          0.231691047482839  % 30
          0.2665056338554453 % 31
          0.3060144741109534 % 32
          0.3508128223486922 % 33
          0.4015694373296012 % 34
          0.4590355489842755 % 35
          0.5240549123983681 % 36
          0.5975750808014414 % 37
          0.6806600449716177 % 38
          0.7745044042722101 % 39
          0.8804492544892654 % 40
          1.0                % 41
    ];
end

for k = 1:1:NLEV_full
    A_full(k,1) = 0.5*( A_half(k,1)+A_half(k+1,1) );
    B_full(k,1) = 0.5*( B_half(k,1)+B_half(k+1,1) );
end
eta_half = A_half + B_half;
eta_full = A_full + B_full;

% 计算地表高度
zs_half = zeros(1, NX_half);
zs_full = zeros(1, NX_full);
right = zs_half(1,2:NX_half);
left  = zs_half(1,1:NX_half-1);
dzsdx_full = (right-left) / delta_X;
        
% 计算地面静力气压
phs_full = ones(1,NX_full) * p_ref; % 没有地形，地面静力气压就是参考气压
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
G = g^2/N_freq^2/Cp;
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        ph_ratio = ph_cell(k,i) / phs_full(i);
        fenzi = T_surf * power(p_ref/phs_full(i), kappa);
        fenmu = T_surf / G * (power(ph_ratio, kappa) - 1) + 1;
        pt(k,i) = fenzi / fenmu;
    end
end

% 计算地势
for k = 1:1:NLEV_half
    for i = 1:1:NX_full
        num_log = T_surf/G* (power(ph_ns(k,i)/phs_full(i), kappa) - 1) + 1;
        z(k,i) = -g/N_freq^2 * log(num_log);
    end
end
gz = g * z;

% determine where to truncate
z_top_set = 10.0e3; % 好像设太大了没法得到？
trunc_lev = 1;
for k = 1:1:NLEV_half
    if z(k,1) < z_top_set
        trunc_lev = k;
        break
    end
end

z_top_calc = z(1,1) %#ok<NOPTS>
ph_top_calc = ph_ns(1,1) %#ok<NOPTS>

% 加上位温扰动
% perturbed pt profile
% 扰动量 = Δpt * sin(pi*z/H)/(1+(x-xc)^2/a^2)
delta_pt = 1; % K
xc = (x_lo+x_hi)/2.0 - u_ref * 3000; 
a = 5.0e3; H = 10.0e3;
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        local_z = (z(k+1,i) + z(k,i)) / 2.0;
        fenzi = sin(pi * local_z / z(1,i));
        fenmu = 1 + power((x_full(i)-xc)/a, 2.0);
        pt_perturb(k,i) = delta_pt * fenzi / fenmu;
    end
end
pt_tot = pt + pt_perturb;

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
simu_time = 1;% 小时
delta_t = 2.0; % 秒
delta_t_to_write = 5; % 分钟
nstep = simu_time * 3600 / delta_t; 
nstep_to_write = delta_t_to_write * 60 / delta_t;
max_delta_t = min_delta_x / (340+u_ref)

kexi = 0.5; % 0.5 for Crank-Nicholson
zh = z_top_calc * 2; % non-hydrostatic
%tau_0 = 1.0 / 0.15; % for dekta_x = 500m 
tau_0 = 1.0 / 0.3;

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

% 输出初始场文件
save("../cmake-build-release/A_half.txt","A_half", "-ascii");
save("../cmake-build-release/B_half.txt","B_half", "-ascii");
save("../cmake-build-release/zs_half.txt","zs_half", "-ascii");
save("../cmake-build-release/zs_full.txt","zs_full", "-ascii");
save("../cmake-build-release/phs0.txt","phs_full", "-ascii");
save("../cmake-build-release/u0.txt","u", "-ascii");
save("../cmake-build-release/w0.txt","w", "-ascii");
save("../cmake-build-release/pt0.txt","pt_tot", "-ascii");
save("../cmake-build-release/geo0.txt","gz", "-ascii");

subplot(2,2,1);
plot(x_full, gz); xlabel("x/m"); ylabel("gz/(m^2/s^2)");
title("geo potential"); 
subplot(2,2,2);
plot(x_full, ph_cell); xlabel("x/m"); ylabel("ph/Pa"); set(gca,'YDir','reverse');
title("hydrostatic pressure");
subplot(2,2,3);
plot(x_full, pt); xlabel("x/m"); ylabel("pt/K");
title("potential temperature");
subplot(2,2,4);
contourf(pt_perturb); xlabel("x/m"); ylabel("perturbed pt/K"); set(gca,'YDir','reverse');
title("perturbed potential temperature");