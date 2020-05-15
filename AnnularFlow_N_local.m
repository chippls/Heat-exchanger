function [Nu, f] =  AnnularFlow_N_local(Re, Pr, ratio, p, relRough)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
delta = 0.01 * ratio;
[Nu_high, f_high] = AnnularFlow_N(Re, Pr, 1.01*ratio, p, relRough);
[Nu_low, f_low] = AnnularFlow_N(Re, Pr, 0.99*ratio, p, relRough);
Nu = (Nu_high * 1.01 - Nu_low * 0.99)*50;    
f = (f_high * 1.01 - f_low * 0.99)*50;
end

