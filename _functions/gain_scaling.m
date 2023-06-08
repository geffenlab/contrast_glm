% gain scaling function
function [lambda_scaled, this_gain] = gain_scaling(obj,lambda,contrast,ops)

this_gain = gain(obj,contrast,ops);
lambda_scaled = log(obj.base_rate) + this_gain .* (lambda-log(obj.base_rate));