function fits = fit2stepGLM(y,X,C,ops)

% clear warnings
lastwarn('');

% if no include specified, use all data points
if ~isfield(ops,'include') || isempty(ops.include)
    ops.include = true(size(y));
end
ops.include = logical(ops.include);

% apply weights (no weighting for now)
ops.weights = ones(size(y));

% apply inclusion
y = y(ops.include);
X = X(ops.include,:);
C = C(ops.include,:);
weights = ops.weights(ops.include);


%% step 1: fit the STRF
tic;
if ops.alpha >= 0 && ops.alpha <= 1
    fprintf('\tStep 1: fitting regularized STRF using glmnet... ');
    
    % use regularization
    opts = glmnetSet;
    opts.alpha = ops.alpha;
    nfolds = 10;
    fit_nogain = cvglmnet(X,y,'poisson',opts,[],nfolds,[],false);
    coeffs_nogain = cvglmnetPredict(fit_nogain,[],'lambda_1se','coefficients');
    
else
    fprintf('\tStep 1: fitting unregularized STRF using glmfit... ');
    
    % initialize strf with OLS estimate
    %if isfield(ops,'initializeSTRF') && ops.initializeSTRF
    %    fprintf(' initializing with OLS estimate... ');
    %    b0 = X \ y;
    %    fit_nogain = fitglm(X,y,'distribution','poisson',...
    %       'weights',weights,...
    %        'B0',[1; b0]);
    %else
    %   fit_nogain = fitglm(X,y,'distribution','poisson',...
    %        'weights',weights);
    %end
    %coeffs_nogain = fit_nogain.Coefficients.Estimate;
        
    [coeffs_nogain, fit_nogain.dev, fit_nogain.stats] = ...
        glmfit(X,y,'poisson');
    
end
toc;

% catch the last warnings
fit_nogain.warning = lastwarn;

% compute strf prediction
preds_nogain = X * coeffs_nogain(2:end); % no intercept!


%% step 2: fit the gain
fprintf('\tStep 2: fit gain coefficients... '); tic;

% scale predictors if specified
if isfield(ops,'scale') && ops.scale
    fprintf(' scaling predictors... ');
    scale = std(preds_nogain);
else
    scale = 1;
end
preds_nogain_scaled = preds_nogain / scale;
dm_gain = [preds_nogain_scaled, preds_nogain_scaled .* C, C];

% fit
[coeffs_gain, fit_gain.dev, fit_gain.stats] = ...
    glmfit(dm_gain,y,'poisson','weights',weights);
fit_gain.pred = glmval(coeffs_gain,dm_gain,'log');

% catch the last warning
fit_gain.warning = lastwarn;

% rescale
coeffs_gain(2:2+size(C,2)-1) = coeffs_gain(2:2+size(C,2)-1) / scale;
toc;


%% format output
fit_nogain.pred = preds_nogain;
fit_nogain.coeffs = coeffs_nogain;
fit_gain.coeffs = coeffs_gain;
fit_gain.scale = scale;

fits.y = y;
fits.strf_fit = fit_nogain;
fits.gain_fit = fit_gain;
    
    
    
    
    
    
    
    
    
    
    