% firing rate function
function [l,linear_drive,g,x0] = lambda(obj,x,contrast,ops)
% x = stimulus

% set the ITI periods where stim == 0 to the stimulus mean, which will
% allow for changes in operating point relative to the centered stimulus
x(x==0) = obj.x0;

% convolve with the strf
x0 = convSTRF(x'-obj.operating_point,fliplr(obj.beta))';

% linear drive
linear_drive = log(obj.base_rate) + x0;

% gain works by rescaling beta. At each lag, the beta at that lag for all
% frequencies is rescaled depending on the contrast at that lag
%g = gain(obj,contrast);
%drive_with_gain = log(obj.base_rate) + ...
%    convSTRF(x'-obj.operating_point,fliplr(obj.beta),'full',g')';
[drive_with_gain,g] = gain_scaling(obj, linear_drive, contrast, ops);

% exponentiate for poisson rate
l = exp(drive_with_gain);
%l(l>200) = 200;
