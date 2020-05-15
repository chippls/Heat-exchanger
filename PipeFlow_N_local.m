function [Nu,f] = PipeFlow_N_local(Re, Pr, ratio, relRough)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if ratio < 0.1
    ratio = 0.1
end
delta = ratio * 0.01;
[Nu_high, f_high] = PipeFlow_N(Re, Pr, 1.01*ratio, relRough);
[Nu_low, f_low] = PipeFlow_N(Re, Pr, 0.99*ratio, relRough);
Nu = (Nu_high * 1.01 - Nu_low * 0.99)*50;    
f = (f_high * 1.01 - f_low * 0.99)*50;
end

