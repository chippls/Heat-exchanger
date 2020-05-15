function [Nu, f] = AnnularFlow_Laminar(Re, Pr, ratio, p)
%ANNULARFLOW_LAMINAR 此处显示有关此函数的摘要
%   此处显示详细说明
    p_m=((1-p^2)/(2*log(1/p)))^(1/2);
    f_fd=(4/Re)*(16*(1-p)^2)/(1+p^2-2*p_m^2);
    x_plus = ratio/Re;
    fR=3.44/sqrt(x_plus)+(1.25/(4*x_plus)+f_fd*Re/4-3.44/sqrt(x_plus))/(1+0.00021*x_plus^(-2));
    f=fR*4/Re;
    x_star=ratio/(Re*Pr);
    Nusselt_T_fd=0.580342564/p+6.09483719-4.45569753*p+2.64812415*p^2;
    DNusselt_T=1.75450933*exp(-0.402783707*log(x_star))^1.05038051;
    if Pr > 0.72
        DNurat=0.6847+0.3153*exp(-1.26544559*(log(Pr)-log(0.72)));
    else
        DNurat=1.68-0.68*exp(0.32*(log(Pr)-log(0.72)));
    end
    Nu = Nusselt_T_fd+DNurat*DNusselt_T;
end

