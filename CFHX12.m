 function [T_h, T_c, P_h, P_c, Nu_h, Nu_c, Re_h, Re_c,T_wo,T_w,Q_h,Q_c,T_c_out,T_h_out,epsion]  = CFHX12(P_h_in,P_c_in,T_h_in, T_c_in, d1,d2,d3,d4, D_coil,cycles, m_h,m_c, roughness,epsilon)
%DIMESIONLESSCFHX 逆流套管式换热器计算,考虑压降，根据CFHX2改造四计算域模型，计算与EES吻合
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
R = 'hydrogen';
m_h = (m_h+m_c)/2;
m_c = m_h;
theta = 5.67*10^(-8);

%计算总管长 m
length = pi * cycles * D_coil;

%划分微元，N为微元数
N = 40;
deltaL = length / N;
T_h = zeros(1,N+1);
T_c = zeros(1,N+1);
T_w = zeros(1,N+2);
T_wo = zeros(1,N+2);

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

coeff_1 = zeros(1,N);
coeff_2 = zeros(1,N);
coeff_3 = zeros(1,N);

Re_h = zeros(1,N);
Re_c = zeros(1,N);
Pr_h = zeros(1,N);
Pr_c = zeros(1,N);
Nu_h = zeros(1,N);
Nu_c = zeros(1,N);
% De_h = zeros(1,N);
h_correct_h = zeros(1,N);
h_correct_c = zeros(1,N);

%计算内外管相对粗糙度
relRough_h = roughness / d1;
relRough_c = roughness / (d3 - d2);

%计算内外流道面积 m^2
A_h = pi * d1^2 / 4;
A_c = pi * (d3^2 - d2^2) / 4;

%计算内外流道单个微元对流换热面积 m^2
deltaA_h = pi * d1 * deltaL;
deltaA_c = pi * d2 * deltaL;
deltaA_co = pi * d3 * deltaL;

%计算内外管轴向导热面积 m^2
A_con = pi * (d2^2 - d1^2) / 4;
A_con_o = pi * (d4^2 - d3^2) / 4;

%计算外管辐射面积
A_rad = pi * d4 * deltaL;
T_a = T_h_in;

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
    T_w(i+1) = T_h_in -i*deltaT/(N+1)+1;
    x_h(i) = (i-0.5) * deltaL;
    x_c(i) = (N+0.5-i) * deltaL;
end

for i = 1:N+2
    T_wo(i) = T_w(i)+1;
end

% eps = 0.01;
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
        cp_h(i) = refpropm('C','T',T_ave_h,'P',P_h(i),R);
        cp_c(i) = refpropm('C','T',T_ave_c,'P',P_c(i),R);
        u_h = m_h / A_h / rho_h;
        u_c = m_c / A_c / rho_c;
        Re_h(i) = rho_h * u_h * D_h_h / mu_h;
        Re_c(i) = rho_c * u_c * D_h_c / mu_c;
        Pr_h(i) = mu_h * cp_h(i) / k_h;
        Pr_c(i) = mu_c * cp_c(i) / k_c;
         [Nu_h(i), ~] = PipeFlow_N_local(Re_h(i), Pr_h(i),x_h(i)/D_h_h, relRough_h);
%         De_h(i) = Re_h(i)*sqrt(d1/D_coil);
%         f_h(i) = 0.3164*Re_h(i)^(-0.25)+0.03*sqrt(d1/D_coil);
%         Nu_h(i) = Pr_h(i)*Re_h(i)*f_h(i)/(8*(1+12.7*sqrt(f_h(i)/8)*(Pr_h(i)^(2/3)-1)));
        [Nu_c(i),~] = AnnularFlow_N_local(Re_c(i), Pr_c(i),x_c(i)/D_h_c, d2/d3, relRough_c);
        
        [~,f_h(i)] = PipeFlow_N_local(Re_h(i), Pr_h(i),x_h(i)/D_h_h, relRough_h);
        [~,f_c(i)] = AnnularFlow_N_local(Re_c(i), Pr_c(i), x_c(i)/D_h_c, d2/d3, relRough_c);
        deltaP_h(i) = f_h(i) *deltaL/(2*D_h_h)*rho_h*u_h^2;
        deltaP_c(i) = f_c(i) *deltaL/(2*D_h_c)*rho_c*u_c^2;
        deltaP_curve_h(i) = deltaP_h(i) * (1+0.0823*(1+d1/D_coil)*(d1/D_coil)^0.53*Re_h(i)^0.25)/1000;
        deltaP_curve_c(i) = deltaP_c(i) * (1+0.0823*(1+(d3-d2)/D_coil)*((d3-d2)/D_coil)^0.53*Re_c(i)^0.25)/1000;
        
        h_h = Nu_h(i) * k_h / D_h_h;
        h_c = Nu_c(i) * k_c / D_h_c;
        h_correct_h(i) = (1+3.6*(1-D_h_h/D_coil)*(D_h_h/D_coil)^0.8)*h_h;
        h_correct_c(i) = (1+3.6*(1-D_h_c/D_coil)*(D_h_c/D_coil)^0.8)*h_c;
        
        k_ss_pre = interp1(M(:,1), M(:,2), (T_w(i)+T_w(i+1))/2, 'liner' );
        k_ss_back = interp1(M(:,1), M(:,2), (T_w(i+2)+T_w(i+1))/2, 'liner' );
        
        k_pre = interp1(M(:,1), M(:,2), (T_wo(i)+T_wo(i+1))/2, 'liner' );
        k_back = interp1(M(:,1), M(:,2), (T_wo(i+2)+T_wo(i+1))/2, 'liner' );
        
        
        a(i) = h_correct_h(i) * deltaA_h;
        b(i) = h_correct_c(i) * deltaA_c;
        c(i) = m_h*cp_h(i)-h_correct_h(i) * deltaA_h/2;
        d(i) = m_c*cp_c(i)-h_correct_c(i) * deltaA_c/2;
        e(i) = h_correct_h(i) * deltaA_h/2 + m_h*cp_h(i);
        f(i) = h_correct_c(i) * deltaA_c/2 + m_c*cp_c(i);


        g(i) = k_ss_pre * A_con / deltaL;
        h(i) = k_ss_back * A_con / deltaL;
        
        %考虑换热器外壁面
        coeff_1(i) = h_correct_c(i) * deltaA_co;
        coeff_2(i) = k_pre*A_con_o / deltaL;
        coeff_3(i) = k_back*A_con_o / deltaL;
        coeff_4 = epsilon*A_rad*theta;
        
        
    end
    T_h_old = T_h;
    T_c_old = T_c;
    T_w_old = T_w;
    T_wo_old = T_wo;

    for i=1:N
        T_h(i+1) = (a(i) * T_w(i+1) + c(i)*T_h(i))/e(i);
        T_c(i) = (b(i)*T_w(i+1)+coeff_1(i)*T_wo(i+1)+(d(i)-coeff_1(i)/2)*T_c(i+1))/(f(i)+coeff_1(i)/2);
%         T_c(i) = (b(i) * T_w(i+1) + d(i)*T_c(i+1))/f(i);
        P_h(i+1) = P_h(i) - deltaP_curve_h(i);
%         P_c(N+1) = P_h(N+1);
        P_c(i+1) = P_c(i) - deltaP_curve_c(i);
        
        if i==1
           T_w(i+1) = (2*g(i)*T_w(i) + h(i)*T_w(i+2) + a(i)*(T_h(i)+T_h(i+1))/2 + b(i)*(T_c(i)+T_c(i+1))/2) / (2*g(i)+h(i)+a(i)+b(i));
           p = [coeff_4,0,0,coeff_3(i)+2*coeff_2(i)+coeff_1(i),-(2*coeff_2(i)*T_wo(i)+coeff_3(i)*T_wo(i+2)+coeff_1(i)*(T_c(i)+T_c(i+1))/2+coeff_4*T_a^4)];
           roots(p);
           %T_wo(i+1) = ans(4);
           T_wo(i+1) = ans(1);
        elseif i==N
           T_w(i+1) = (g(i)*T_w(i) + 2*h(i)*T_w(i+2) + a(i)*(T_h(i)+T_h(i+1))/2 + b(i)*(T_c(i)+T_c(i+1))/2) / (g(i)+2*h(i)+a(i)+b(i));
           p = [coeff_4,0,0,2*coeff_3(i)+coeff_2(i)+coeff_1(i),-(coeff_2(i)*T_wo(i)+2*coeff_3(i)*T_wo(i+2)+coeff_1(i)*(T_c(i)+T_c(i+1))/2+coeff_4*T_a^4)];
           roots(p);
           %T_wo(i+1) = ans(4);
           T_wo(i+1) = ans(1);
        else
           T_w(i+1) = (g(i)*T_w(i) + h(i)*T_w(i+2) + a(i)*(T_h(i)+T_h(i+1))/2 + b(i)*(T_c(i)+T_c(i+1))/2) / (g(i)+h(i)+a(i)+b(i));
           p = [coeff_4,0,0,coeff_3(i)+coeff_2(i)+coeff_1(i),-(coeff_2(i)*T_wo(i)+coeff_3(i)*T_wo(i+2)+coeff_1(i)*(T_c(i)+T_c(i+1))/2+coeff_4*T_a^4)];
           roots(p);
           %T_wo(i+1) = ans(4);
           T_wo(i+1) = ans(1);
        end
       
        
        
    end
%     T_wo(1) = T_w(2);
%     T_wo(N+2) = T_wo(N+1);
    
T_wo(1) = T_c(1);
T_wo(N+2) = T_c(N+1);
% T_w(1) = (T_h(1)+T_c(1))/2;
% T_w(N+2) = (T_h(N+1)+T_c(N+1))/2;
% T_w(1) = T_h(1);
% T_w(N+2) = T_h(N+1);


    j = j+1
    T_old = [T_h_old,T_c_old,T_w_old];
    T = [T_h,T_c,T_w];
    if(stopCiteria(T,T_old))
        break;
    end
    
    h_enthalpy_in = refpropm('H','T',T_h(1),'P',P_h(1),R);
    h_enthalpy_out = refpropm('H','T',T_h(N+1),'P',P_h(N+1),R);
    c_enthalpy_in = refpropm('H','T',T_c(N+1),'P',P_c(N+1),R);
    c_enthalpy_out = refpropm('H','T',T_c(1),'P',P_c(1),R);
    h_ideal_out = refpropm('H','T',T_c(N+1),'P',P_h(N+1),R);
    c_ideal_out = refpropm('H','T',T_h(1),'P',P_c(1),R);
    Q_h = m_h*(h_enthalpy_in-h_enthalpy_out);
    Q_c = m_c*(c_enthalpy_out-c_enthalpy_in);
    Q_h_ideal = m_h*(h_enthalpy_in-h_ideal_out);
    Q_c_ideal = m_c*(c_ideal_out-c_enthalpy_in);
    epsion = min(Q_h,Q_c)/min(Q_h_ideal,Q_c_ideal)  
    T_h_out = T_h(N+1)
    T_c_out = T_c(1)

end

%考虑涡流影响的网格尺寸
eta_h = u_h^3/length;
eta_c = u_c^3/length;
% yeta_h = ((mu_h/rho_h)^3/eta_h)^0.25
% yeta_c = ((mu_c/rho_c)^3/eta_c)^1/4

% dx = linspace(0,1,N+1);
% plot(dx,T_h,dx,T_c);
% T_h_out = T_h(N+1)
% T_c_out = T_c(1)


