function [T_h, T_c, T_w, Nu_h, Nu_c, Re_h, Re_c,P_h, P_c] = CFHX3(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_dot, roughness)
%DIMESIONLESSCFHX 逆流套管式换热器计算
%   此处显示详细说明
%   P_h_in      换热器高压侧入口压力       kPa
%   P_c_in      换热器低压侧入口压力       kPa
%   T_h_in      换热器高压侧入口温度       K
%   T_c_in      换热器低压侧入口温度       K
%   d1          换热器内管内径             m
%   d2          换热器内管外径             m
%   d3          换热器外管内径             m
%   d4          换热器外管外径             m
%   D_coil      螺旋管大径                 m
%   cycles      螺旋管盘绕圈数 
%   m_dot       质量流量                   kg/s 
%   roughness   管子粗糙度                 m




M = load('SS304.txt');
R = 'helium';
m_h = m_dot;
m_c = m_dot;

%计算总管长 m
length = pi * cycles * D_coil;

%划分微元，N为微元数
N = 100;
deltaL = length / N;
T_h = zeros(1,N+1);
T_c = zeros(1,N+1);
T_w = zeros(1,N+2);

P_h = zeros(1,N+1);
P_c = zeros(1,N+1);
f_h = zeros(1,N);
f_c = zeros(1,N);
deltaP_h = zeros(1,N);
deltaP_c = zeros(1,N);
deltaP_curve_h = zeros(1,N);
deltaP_curve_c = zeros(1,N);


x_h = zeros(1,N);
x_c = zeros(1,N);
a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,N);
d = zeros(1,N);
e = zeros(1,N);
f = zeros(1,N);
g = zeros(1,N);
h = zeros(1,N);

a_o = zeros(1,N);
b_o = zeros(1,N);
c_o= zeros(1,N);
d_o= zeros(1,N);
e_o= zeros(1,N);
T_wo = zeros(1,N);
Q_rad = zeros(1,N);

Re_h = zeros(1,N);
Re_c = zeros(1,N);
Pr_h = zeros(1,N);
Pr_c = zeros(1,N);
Nu_h = zeros(1,N);
Nu_c = zeros(1,N);
% De_h = zeros(1,N);

%计算内外管相对粗糙度
relRough_h = roughness / d1;
relRough_c = roughness / (d4 - d3);

%计算内外流道面积 m^2
A_h = pi * d1^2 / 4;
A_c = pi * (d3^2 - d2^2) / 4;

%计算内外流道单个微元对流换热面积 m^2
deltaA_h = pi * d1 * deltaL;
deltaA_c = pi * d2 * deltaL;
deltaA_co = deltaA_c*(d3/d2);

%计算内外管轴向导热面积 m^2
A_con = pi * (d2^2 - d1^2) / 4;
A_con_o = pi * (d4^2 - d3^2) / 4;

%外壁面辐射漏热部分
tao = 0.18;
garma = 5.67e-8;
T_amb = 300;
A_rad = pi*d4*deltaL;

%计算内外流道水利直径 m
D_h_h = d1;
D_h_c = d3 - d2;

%计算最大温差 K
deltaT = T_h_in - T_c_in;

T_h(1) = T_h_in;
T_c(N+1) = T_c_in;

%初始化温度分布 线性分布
for i = 1:N+1
    T_h(i) = T_h_in - (i-1)*deltaT / N;
    T_c(i) = T_h(i);
    
    P_h(i) = P_h_in;
    P_c(i) = P_c_in;
    
end
T_w(1) = T_h_in;
T_w(N+2) = T_c_in;
for i = 1:N
    T_w(i+1) = T_h_in -i*deltaT/(N+1);
    x_h(i) = (i-0.5) * deltaL;
    x_c(i) = (N+0.5-i) * deltaL;
end
for i = 1:N+2
    T_wo(i) = T_w(i);
end
eps = 0.01;
j = 1;

while 1==1
    for i = 1:N
        T_ave_h = (T_h(i)+T_h(i+1)) / 2;
        T_ave_c = (T_c(i)+T_c(i+1)) / 2;
        rho_h = refpropm('D','T',T_ave_h,'P',P_h(i),R);
        rho_c = refpropm('D','T',T_ave_c,'P',P_c(i),R);
        k_h = refpropm('L','T',T_ave_h,'P',P_h(i),R);
        k_c = refpropm('L','T',T_ave_c,'P',P_c(i),R);
        mu_h = refpropm('V','T',T_ave_h,'P',P_h(i),R);
        mu_c = refpropm('V','T',T_ave_c,'P',P_c(i),R); 
        cp_h = refpropm('C','T',T_ave_h,'P',P_h(i),R);
        cp_c = refpropm('C','T',T_ave_c,'P',P_c(i),R);
        u_h = m_h / A_h / rho_h;
        u_c = m_c / A_c / rho_c;
        Re_h(i) = rho_h * u_h * D_h_h / mu_h;
        Re_c(i) = rho_c * u_c * D_h_c / mu_c;
        Pr_h(i) = mu_h * cp_h / k_h;
        Pr_c(i) = mu_c * cp_c / k_c;
        [Nu_h(i), ~] = PipeFlow_N_local(Re_h(i), Pr_h(i),x_h(i)/D_h_h, relRough_h);
%         De_h(i) = Re_h(i)*sqrt(d1/D_coil);
%         Nu_h(i) = 0.0509*Re_h(i)^0.817*Pr_h(i)^(0.3)*(d1/D_coil)^(-0.1);
        [Nu_c(i),~] = AnnularFlow_N_local(Re_c(i), Pr_c(i),x_c(i)/D_h_c, d2/d3, relRough_c);
        
        [~,f_h(i)] = PipeFlow_N_local(Re_h(i), Pr_h(i),x_h(i)/D_h_h, relRough_h);
        [~,f_c(i)] = AnnularFlow_N_local(Re_c(i), Pr_c(i), x_c(i)/D_h_c, d2/d3, relRough_c);
        deltaP_h(i) = f_h(i) *deltaL/(2*D_h_h)*rho_h*u_h^2;
        deltaP_c(i) = f_c(i) *deltaL/(2*D_h_c)*rho_c*u_c^2;
        deltaP_curve_h(i) = deltaP_h(i) * (1+0.0823*(1+d1/D_coil)*(d1/D_coil)^0.53*Re_h(i)^0.25)/1000;
        deltaP_curve_c(i) = deltaP_c(i) * (1+0.0823*(1+(d3-d2)/D_coil)*((d3-d2)/D_coil)^0.53*Re_c(i)^0.25)/1000;
        
        h_h = Nu_h(i) * k_h / D_h_h;
        h_c = Nu_c(i) * k_c / D_h_c;
        h_correct_h = (1+3.6*(1-D_h_h/D_coil)*(D_h_h/D_coil)^0.8)*h_h;
        h_correct_c = (1+3.6*(1-D_h_c/D_coil)*(D_h_c/D_coil)^0.8)*h_c;
        
        k_ss_pre = interp1(M(:,1), M(:,2), (T_w(i)+T_w(i+1))/2, 'liner' );
        k_ss_back = interp1(M(:,1), M(:,2), (T_w(i+2)+T_w(i+1))/2, 'liner' );
        
        k_po = interp1(M(:,1), M(:,2), (T_wo(i)+T_wo(i+1))/2, 'liner' );
        k_bo = interp1(M(:,1), M(:,2), (T_wo(i+2)+T_wo(i+1))/2, 'liner' );
        
        a(i) = h_correct_h * deltaA_h;
        b(i) = h_correct_c * deltaA_c;
        c(i) = m_h*cp_h-h_correct_h * deltaA_h/2;
%         d(i) = m_c*cp_c-h_correct_c * deltaA_c/2;
        d(i) = m_c*cp_c-h_correct_c * (deltaA_c+deltaA_co)/2;
        e(i) = h_correct_h * deltaA_h/2 + m_h*cp_h;
%         f(i) = h_correct_c * deltaA_c/2 + m_c*cp_c;
        f(i) = h_correct_c * (deltaA_c+deltaA_co)/2 + m_c*cp_c;
        b_o(i) = h_correct_c*deltaA_co;
        

        g(i) = k_ss_pre * A_con / deltaL;
        h(i) = k_ss_back * A_con / deltaL;
        a_o(i) = (k_po+k_bo)*A_con_o/deltaL+h_correct_c*deltaA_co;
        c_o(i) =  k_po*A_con_o/deltaL;
        d_o(i) = k_bo*A_con_o/deltaL;
        e_o(i) = h_correct_c*deltaA_co;
        
        
    end
    T_h_old = T_h;
    T_c_old = T_c;
    T_w_old = T_w;
    T_wo_old = T_wo;
    for i=1:N
        T_h(i+1) = (a(i) * T_w(i+1) + c(i)*T_h(i))/e(i);
        T_c(i) = (b(i) * T_w(i+1) + d(i)*T_c(i+1)+b_o(i)*T_wo(i+1))/f(i);
        P_h(i+1) = P_h(i) - deltaP_curve_h(i);
        P_c(N+1) = P_h(N+1);
        P_c(i) = P_c(i+1) - deltaP_curve_c(i);
        Q_rad(i+1) = tao*garma*(T_amb^4-T_wo(i+1)^4)*A_rad;
        
        if i==1
           T_w(i+1) = (2*g(i)*T_w(i) + h(i)*T_w(i+2) + a(i)*(T_h(i)+T_h(i+1))/2 + b(i)*(T_c(i)+T_c(i+1))/2) / (2*g(i)+h(i)+a(i)+b(i));
           T_wo(i+1) = (2*c_o(i)*T_wo(i)+d_o(i)*T_wo(i+2)+e_o(i)*(T_c(i)+T_c(i+1))/2+Q_rad(i+1))/(a_o(i)+c_o(i));
        elseif i==N
           T_w(i+1) = (g(i)*T_w(i) + 2*h(i)*T_w(i+2) + a(i)*(T_h(i)+T_h(i+1))/2 + b(i)*(T_c(i)+T_c(i+1))/2) / (g(i)+2*h(i)+a(i)+b(i));
           T_wo(i+1) = (c_o(i)*T_wo(i)+d_o(i)*T_wo(i+2)+e_o(i)*(T_c(i)+T_c(i+1))/2+Q_rad(i+1))/a_o(i);
        else
           T_w(i+1) = (g(i)*T_w(i) + h(i)*T_w(i+2) + a(i)*(T_h(i)+T_h(i+1))/2 + b(i)*(T_c(i)+T_c(i+1))/2) / (g(i)+h(i)+a(i)+b(i));
           T_wo(i+1) = (c_o(i)*T_wo(i)+2*d_o(i)*T_wo(i+2)+e_o(i)*(T_c(i)+T_c(i+1))/2+Q_rad(i+1))/(a_o(i)+d_o(i));
        end
        
        
    end
    T_w(1) = T_w(2);
    T_w(N+2) = T_w(N+1);
    T_wo(1) = T_wo(2);
    T_wo(N+2) = T_wo(N+1);
%     T_w(1) = (T_h(1)+T_c(1))/2
%     T_w(N+2) = (T_h(N+1)+T_c(N+1))/2
    T_h_error = abs(T_h-T_h_old);
    T_c_error = abs(T_c-T_c_old);
    T_w_error = abs(T_w-T_w_old);
    T_wo_error = abs(T_wo-T_wo_old);
    
    max_th = max(T_h_error);
    max_tc = max(T_c_error);
    max_tw = max(T_w_error);
    max_two = max(T_wo_error);
    M_error = [max_th,max_tc,max_tw,max_two];
    max_error = max(M_error)
    j = j+1
    if(max_error < eps)
        break;
    end
    
%     T_h = T_h_old/2 + T_h/2;
%     T_c = T_c_old/2 + T_c/2;
%     T_w = T_w_old/2 + T_w/2;
    
end


