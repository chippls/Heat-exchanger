% P_h_in = 0.4*1000;
% P_c_in = 0.1*1000;
% T_h_in = 305;
% T_c_in = 74;
% d1 = 1.753/1000;
% d2 = 3.175/1000;
% d3 = 4.752/1000;
% d4 = 6.35/1000;
% D_coil = 30/1000;
% cycles = 9;
% m_dot = 5.6/1000000;
% roughness = 0.000001;
clc

P_h_in = 0.405*1000;
P_c_in = 0.404*1000;
T_h_in = 282.65;
T_c_in = 89.46;
d1 = 1.753/1000;
d2 = 3.175/1000;
d3 = 4.752/1000;
d4 = 6.35/1000;
D_coil = 30/1000;
cycles = 8;
m_dot = 5.83/1000000;
roughness = 0.000001;

% [T_h, T_c, T_w] = dimesionlessCFHX(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_dot, roughness)
[T_h, T_c, T_w, Nu_h, Nu_c, Re_h, Re_c,P_h,P_c] = CFHX12(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_dot, roughness);
result_1 = [ Nu_h; Nu_c; Re_h; Re_c];
result_2 = [T_h;T_c;P_h;P_c];
xlswrite('result.xlsx',result_1);
xlswrite('result.xlsx',result_2,2);
