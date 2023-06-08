function [res,ops] = fitGLM(X,c,y,ops)
    
%% function [res,ops] = fitGLM(X,C,y,ops)
%
% fits a GLM to data to estimate varying gain underlying the
% contrast of a stimulus
%
% x - stimulus design matrix
% c - contrast vector
% y - spikes
% ops - struct for fitting options, relevant options below
%
% res - output struct of results
%
% ops.include - logical vector of timepoints to include in the model
% ops.weights - vector of weights for each timepoint (only valid
%               for non-regularized fits using glmfit)  
% ops.scale - if true, z-scores stimulus predictors
% ops.alpha - single value or vector that specifies regularization
%             alpha<0 = no regularization, alpha>0&alpha<1 = regularized
% ops.dummycode - value for splitting out predictors:
%                 0 = no split, 1 = split only b2, 2 = split b2 & b3
    
tic;  
    
%% SETUP
% contrast predictors
cc = ops.mean_contrast ./ c;
C = lagDesignMatrix(cc,ops.gain_lags)';

% index the predictors and spike rate
X = X(ops.include,:);
C = C(ops.include,:);
y = y(ops.include);

% scale predictors
ops.scale = false;
if ops.scale
    fprintf('Z-scored predictors...\n');
end

% enforce weights
if ~isfield(ops,'weights')

    % without weights, set all to 1
    ops.weights = ones(size(y));
    
end




%% first step: fit the STRF

if ops.alpha(1) >= 0 && ops.alpha(1) <= 1
    fprintf('Step 1: fitting STRF using glmnet... ');
    
    % use regularization
    opts = glmnetSet;
    opts.alpha = ops.alpha(1);
    nfolds = 10;
    fit_nogain = cvglmnet(X,y,'poisson',opts,[],nfolds,[],false);
    coeffs_nogain = cvglmnetPredict(fit_nogain,[],'lambda_1se','coefficients');
    
else
    fprintf('Step 1: fitting STRF using glmfit... ');
    
    % unregularized
    [coeffs_nogain, fit_nogain.dev, fit_nogain.stats] = ...
        glmfit(X,y,'poisson','weights',ops.weights);
    
end
preds_nogain = X * coeffs_nogain(2:end); % no intercept!


toc;


%% second step: fit the gain

tic;

% dummy coded contrast
C1 = C;
if ops.dummycode > 0
    
    % ops.shift = [];
    if ~isempty(ops.shift)
        % shift to cover transition time
        cs = circshift(cc,ops.shift);
    else
        cs = cc;
    end
    
    % split contrast predictors
    uc = unique(cc,'stable');
    c1 = C; c2 = C;
    c1(cs == uc(2),:) = 0;
    c2(cs == uc(1),:) = 0;
    
    if ops.dummycode == 2

        % dummy code beta3
        C1 = [c1 c2];
        
    end
    C = [c1 c2];
    
end

% generate contrast and stimulus predictors
XC = preds_nogain .* C;
%XC(C==1) = 1;
dm_gain = [preds_nogain, XC, C1];

plot_on = false;
if plot_on
    figure;
    
    ax1 = subplot(1,3,1);
    clims = linspace(-max(abs(C(:))),max(abs(C(:))),1000);
    cmap = zeroCMap(clims,0);
    imagesc(C(1:300,:)); 
    h = colorbar; colormap(ax1,cmap);
    caxis([clims(1) clims(end)]);
    title('C');
    
    ax2 = subplot(1,3,2);
    clims = linspace(-max(abs(XC(:))),max(abs(XC(:))),1000);
    cmap = zeroCMap(clims,0);
    imagesc(XC(1:300,:)); 
    h = colorbar; colormap(ax2,cmap);
    caxis([clims(1) clims(end)]);
    title('XC');
    
    ax3 = subplot(1,3,3);
    clims = linspace(-max(abs(dm_gain(:))),max(abs(dm_gain(:))),1000);
    cmap = zeroCMap(clims,0);
    imagesc(dm_gain(1:300,:)); 
    h = colorbar; colormap(ax3,cmap);
    caxis([clims(1) clims(end)]);
    title('dm_gain','interpreter','none');
end

if ops.scale
    scale = std(preds_nogain);
    preds_nogain_scaled = preds_nogain ./ scale;
    dm_gain = [preds_nogain_scaled, preds_nogain_scaled .* C, C1];
end  

% fit
if ops.alpha(2) >= 0 && ops.alpha(2) <= 1
    fprintf('Step 2: fitting gain using glmnet... ');
    
    % use regularization
    opts = glmnetSet;
    opts.alpha = ops.alpha(2);
    nfolds = 10;
    fit_gain = cvglmnet(dm_gain,y,'poisson',opts,[],nfolds,[],false);
    coeffs_gain = cvglmnetPredict(fit_gain,[],'lambda_1se','coefficients');
    fit_gain.pred = cvglmnetPredict(fit_gain,dm_gain,'lambda_1se','response');
    
else
    fprintf('Step 2: fitting gain using glmfit... ');
    
    % unregularized
    [coeffs_gain, fit_gain.dev, fit_gain.stats] = ...
        glmfit(dm_gain,y,'poisson','weights',ops.weights);
    fit_gain.pred = glmval(coeffs_gain,dm_gain,'log');
    fit_gain.corr = corr(fit_gain.pred,y);
    
end

% goodness of fit
fit_gain.corr = corr(fit_gain.pred,y);

if ops.scale
    % rescale
    coeffs_gain(2:2+size(XC,2)-1) = coeffs_gain(2:2+size(XC,2)-1) ./ scale';
end


toc;


%% post processing
fprintf('Post processing data... '); tic;

% second step coefficients
beta0 = coeffs_gain(1);
beta1 = coeffs_gain(2);
beta2 = coeffs_gain(3:3+size(XC,2)-1);
beta3 = coeffs_gain(3+size(XC,2):end);

% old gain index (when C isn't masked)
%w_norm = (beta1+C*beta2) / (beta1+sum(beta2));

% gain index (ratio of gradient of changing contrast model to
% gradient of a model with constant contrast)
C_new = C;
C_new(C ~= 0) = 1; % constant contrast model
w_norm = (beta1+C*beta2) ./ (beta1+C_new*beta2);
w = nan(size(ops.include));
w(ops.include) = w_norm;

% results structure
res.strf_fit = fit_nogain;
res.strf_fit.coeffs = coeffs_nogain;
res.strf_fit.pred = preds_nogain;
res.strf = reshape(res.strf_fit.coeffs(2:end),length(ops.f),[]);
res.gain_fit = fit_gain;
res.beta0 = beta0;
res.beta1 = beta1;
res.beta2 = beta2;
res.beta3 = beta3;
res.w = w;


toc;