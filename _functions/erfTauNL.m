function err = erfTauNL(x,y,tau,res,ops)
% function err = erfTauNL(x,y,tau,res,ops)
%  x: STRF convolved with stimulus, trial mean
%  y: firing rate
%  tau: time constant in samples
%  res: results structure
%  ops: options from results

if ~isfield(ops,'valid_conv')
    warning('erfTauNL.m: ops.valid_conv doesnt exist, using whole trace');
    valid_conv = 1;
else
    valid_conv = ops.valid_conv;
end

% calculate firing rate with adapting NL
yh = tauNL(x,tau,res,ops);

% calculate error where convolution is valid
err = norm(yh(valid_conv:end)-y(valid_conv:end));
