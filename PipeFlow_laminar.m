function [Nu, f] = PipeFlow_laminar(Re, Pr, ratio)
%PIPEFLOW_LAMINAR 此处显示有关此函数的摘要
%   此处显示详细说明
Gz = Re * Pr / ratio;
x = ratio / Re;
fR = 3.44/sqrt(x)+(1.25/(4*x)+16-3.44/sqrt(x))/(1+0.00021*x^(-2));
f=4*fR/Re;
Gm=Gz^(1/3);
Nu = 3.66+((0.049+0.02/Pr)*Gz^1.12)/(1+0.065*Gz^0.7);
end

