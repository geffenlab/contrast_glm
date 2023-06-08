function [y,p_tau] = tauNL(x,tau,res,ops)
% function [y,p_tau] = tauNL(x,tau,res,ops)
%  x: STRF convolved with stimulus, trial mean
%  tau: time constant in samples
%  res: results structure
%  ops: options from results

% normalized contrast
nc = normalize(res.c,'range',[0 1]);

% time vector (max length is the block length)
T = 1:(ops.blockLength * ops.fs);

% exponential kernel
e_kernel = exp(-T/tau);

% smoothly vary parameters according to tau
p = res.gc_ln.ahat;
for i = 1:size(p,2)
    tmp = p(1,i) + ( p(2,i) - p(1,i) ) .* (conv(e_kernel,nc)/sum(e_kernel));
    p_tau(i,:) = tmp(1:length(nc));
end

% generate rate prediction at each timestep as the parameters adapt
mdl = ops.ln_model;
for t = 1:length(x)
    y(t) = mdl(p_tau(:,t),x(t));
end
