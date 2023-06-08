% simulation
function [y,l,h,g,x0] = simulate(obj,x,ops)
% y = spike count
% l = spike rate
% h = linear drive
% g = gain
contrast = obj.contrast;
[l,h,g,x0] = lambda(obj,x,contrast,ops);
y = poissrnd(l);