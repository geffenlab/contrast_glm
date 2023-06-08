function w = computew(beta1,beta2,C,type,ns)

%% function w = computew(beta1,beta2,C,type,ns)
%
% this function computes different versions of w, which is an
% estimate of the timecourse of neural gain computed from
% parameters fit using a GLM that is sensitive to multiplicative
% effects of stimulus variability (contrast)
%
% INPUTS
%  beta1: parameters scaling the stimulus response (TF x 1)
%  beta2: parameters scaling the stimulus x contrast interaction (TF x 1)
%  C: design matrix for contrast (TF x n timepoints)
%  type: selector for different estimates of w, options are:
%  ns: number of samples per trial, to reshape the time course of w (optional)

if strcmp(type,'simple')
    
    % original formula, based on inverting the GLM to estimate
    % parameters of the forward model
    bs = beta1 + beta2;
    bs(bs==0) = 1e-5; % avoid divide by 0
    w = 1 + (beta2'./bs'*(C-1)) ./ numel(beta1);
    
elseif strcmp(type,'covx')
    
    % new version, which utilizes the covariance of x to compare
    % a model with gain control and one without, so is essentially
    % comparing model predictions, not the parameters within the model
    gamma = beta2 .* C;
    dprod = (beta1 + beta2)' * (beta1+gamma);
    bbnorm = norm(beta1+beta2)^2;
    bgnorm = sum((beta1+gamma).^2,1);
    w = 2 * dprod ./ ...
        (bbnorm - bgnorm + ...
         sqrt((bbnorm+bgnorm).^2 - 4*(bbnorm*bgnorm-dprod.^2)));
    
elseif strcmp(type,'norm')
    
    % newest version, a generalized version of the 'simple'
    % version, extended from the original scalar model to 2d model
    gamma_s = beta2 .* C; % TF x n_timepoints
    delta = gamma_s - beta2;
    w = 1 + sum((beta1+beta2).*delta, 1)/(norm(beta1+beta2)^2);
    
end

% if we know the number samples per trial, reshape
if exist('ns','var') & ~isempty(ns)
    
    w = reshape(w,ns,[]')';
    
end    
    
