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

% %氢试验台数据
% P_h_in = 0.982*1000;
% P_c_in = 0.109*1000;
% T_h_in = 299.6;
% T_c_in = 20.89;
% m_h = 30.1/1000000;
% m_c = 29.4/1000000;
d1 = 1.753/1000;
d2 = 3.175/1000;
d3 = 7.747/1000;
d4 = 9.525/1000;
D_coil = 200/1000;
cycles = 20.6;
% clear
% clc

%思科优化后0.2mm壁厚套管
% d1 = 1.6/1000;
% d2 = 2/1000;
% d3 = 3.1/1000;
% d4 = 3.5/1000;
% P_h_in = 0.402*1000;
% P_c_in = 0.40*1000;
% T_h_in = 303.08;
% T_c_in = 96.54;
% cycles = 10;
% m_dot = 

%思科未优化标准管换热器
% P_h_in = 0.405*1000;
% P_c_in = 0.404*1000;
% T_h_in = 282.65;
% T_c_in = 89.46;
% d1 = 1.753/1000;
% d2 = 3.175/1000;
% d3 = 4.752/1000;
% d4 = 6.35/1000;
% D_coil = 30/1000;
% m_dot = 5.83/1000000;
% cycles = 8;
% m_h = m_dot;
% m_c = m_dot;
roughness = 0.000001;
%epsilon = 0.00001;
epsilon = 0;  

% [T_h, T_c, P_h, P_c, Nu_h, Nu_c, Re_h, Re_c]  = CFHX2(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_dot, roughness)
%读取实验数据：高压进口温度K，低压进口温度K，高、低压压力bar，高低压流量mg/s
% data=xlsread("E:/液氢实验/实验数据/开式液氢JT制冷机实验/液氢实验工况汇总.xlsx",'B2:I11');
% NumsOfExp = size(data);
% for i=1:NumsOfExp(1)
%     T_h_in = data(i,1);
%     T_c_in = data(i,3);
%     m_h = data(i,5)/1000000;
%     m_c = data(i,6)/1000000;
%     P_h_in = data(i,7)*101;
%     P_c_in = data(i,8)*101;
for i=1:5
    P_h_in = 8.2*101;
    P_c_in = 101*1.1;
    m_h = (5+5*i)/1000000;
    m_c = m_h;
    T_h_in = 283;
    T_c_in = 21;
    [T_h, T_c, P_h, P_c, Nu_h, Nu_c, Re_h, Re_c,T_wo,T_w,Q_h,Q_c,T_c_out,T_h_out,epsion] = CFHX12(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_h,m_c, roughness, epsilon);
%     [T_h, T_c, P_h, P_c, Nu_h, Nu_c, Re_h, Re_c,T_h_out,T_c_out]  = CFHX2(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_h,m_c, roughness)
    result_1 = [ Nu_h; Nu_c; Re_h; Re_c]';
    result_2 = [T_h_out,T_c_out,epsion];
    result_3 = [T_w;T_wo]';
    start_pos = strcat('K',num2str(i+1));
    end_pos = strcat('M',num2str(i+1));
    xlswrite("E:/液氢实验/实验数据/开式液氢JT制冷机实验/液氢模拟.xlsx",result_2,strcat(start_pos,':',end_pos));
end
%     xlswrite("E:\2020工程热物理数据\无外壁面.xlsx",result_2,1,'A2');

