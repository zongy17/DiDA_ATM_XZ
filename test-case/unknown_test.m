clear
% Unknown test coeffs 

NX_full = 30; NX_half = NX_full + 1;
NLEV_full = 26; NLEV_half = NLEV_full + 1;

% constants
g = 9.80616; 
gamma = 1.4; Cp = 1004.0; Rd = 287.04;
% mount profile: linear mountain, unit: m
% h = hc * ac^2 / (ac^2 + (x-xc)^2)
hc = 800.0; xc = 0.0; ac = 2.0e3;
% ref state
u_ref = 10.0; w_ref = 0.0; 
T_ref = 300.0; p_ref = 1.0e5;
N_freq = g / sqrt(Rd * T_ref);
% domain size
x_lo = -15.0e3; x_hi = 15.0e3;
z_top = 25.0e3;
ph_top = p_ref * exp(- N_freq*N_freq*z_top / g);

% 计算水平坐标
delta_X = (x_hi-x_lo)/NX_full;
x_half = x_lo : delta_X : x_hi;
for i = 1:1:NX_full
    x_full(i) = (x_half(i)+x_half(i+1))/ 2.0;
end

% 计算垂直坐标
A_half =[   0.002194067, 0.004895209, 0.009882418, 0.01805201,...
            0.02983724, 0.04462334 , 0.06160587, 0.07851243,...
            0.07731271, 0.07590131 , 0.07424086, 0.07228744,...
            0.06998933, 0.06728574, 0.06410509, 0.06036322,...
            0.05596111, 0.05078225, 0.04468960, 0.03752191,...
            0.02908949, 0.02084739, 0.01334443, 0.00708499,...
            0.00252136, 0.00000000, 0.00000000];
B_half = [  0.000000, 0.000000, 0.000000, 0.000000,...
            0.000000, 0.000000, 0.000000, 0.000000,...
            0.01505309, 0.03276228, 0.05359622, 0.07810627,...
            0.1069411,  0.1408637, 0.1807720, 0.2277220,...
            0.2829562, 0.3479364, 0.4243822, 0.5143168,...
            0.6201202, 0.7235355, 0.8176768, 0.8962153,...
            0.9534761, 0.9851122, 1.0000000];
eta_half = A_half + B_half;
for k = 1:NLEV_full
    A_full(k) = 0.5 * (A_half(k) + A_half(k+1));
    B_full(k) = 0.5 * (B_half(k) + B_half(k+1));
    eta_full(k) = A_full(k) + B_full(k);
end

% 计算地表高度
zs_half = hc ./ (1 + (x_half-xc)/ac .* (x_half-xc)/ac);
zs_full = hc ./ (1 + (x_full-xc)/ac .* (x_full-xc)/ac);

% 计算地面静力气压
phs_full = p_ref * exp(- g * zs_full / Rd / T_ref);

% 计算各个位置full level的静力气压
for k = 1:1:NLEV_full
    for i = 1:1:NX_full
        ph_cell(k,i) = A_full(k) * p_ref + B_full(k) * phs_full(i);
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
pt = T_ref * power(ph_cell / p_ref, (1-gamma) / gamma);

% 计算地势
gz = g * g / (N_freq*N_freq) * log(p_ref./ph_ns);


subplot(1,2,1);
plot(x_full, gz); xlabel("x/m"); ylabel("gz/(m^2/s^2)");
title("geo potential"); 
subplot(1,2,2);
plot(x_full, ph_cell); xlabel("x/m"); ylabel("ph/Pa");
title("hydrostatic pressure"); 

        
                                