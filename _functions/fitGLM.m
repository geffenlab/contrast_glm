function [res,ops] = fitGLM(X0,y0,cc0,ops)
    
%% function [res,ops] = fitGLM(X0,y0,cc0,ops)
%
% fits a GLM to data to estimate varying gain underlying the
% contrast of a stimulus
%
% X0 - stimulus design matrix
% cc0 - contrast vector
% y0 - spikes
% ops - struct for fitting options, relevant options below
%
% res - output struct of results
%
% ops.include - logical vector of timepoints to include in the model
% ops.alpha - single value or vector that specifies regularization
%             alpha<0 = no regularization, alpha>0&alpha<1 = regularized
% ops.dummycode - value for splitting out predictors:
%                 0 = no split, 1 = split only b2, 2 = split b2 & b3
% ops.cvfolds - number of folds for crossvalidation (set to 1 to skip)
    
tic;  
    
%% SETUP
% original predictors
C0 = lagDesignMatrix(cc0,ops.gain_lags)';

% crossvalidation
cv_trials = crossvalind('kfold',ops.order_r(:,1),ops.cvfolds);
cv_sample = cleanResample(cv_trials',ops.blockLength,ops.period)';

fprintf(sprintf('Fitting 2step GLM with %d folds: ',ops.cvfolds));
tic;
for i = 1:ops.cvfolds
    
    if ops.cvfolds > 1
        trainI = cv_sample ~= i;
        testI = cv_sample == i;
    else
        trainI = cv_sample == i;
        testI = trainI;
    end
    
    % index the predictors and spike rate
    X = X0(ops.include & trainI,:);
    C = C0(ops.include & trainI,:);
    y = y0(ops.include & trainI);
    cc = cc0(ops.include & trainI);
    
    % dummy code the contrast
    [C,C1] = dummycode(C,cc,ops);

    %% first step: fit the STRF

    if ops.alpha(1) >= 0 && ops.alpha(1) <= 1
        %fprintf('Step 1: fitting STRF using glmnet... ');
        
        % use regularization
        opts = glmnetSet;
        opts.alpha = ops.alpha(1);
        nfolds = 10;
        fit_nogain(i) = cvglmnet(X,y,'poisson',opts,[],nfolds,[],false);
        coeffs_nogain(i,:) = cvglmnetPredict(fit_nogain(i),[],'lambda_1se','coefficients');
        
    else
        %fprintf('Step 1: fitting STRF using glmfit... ');
        
        % unregularized
        [coeffs_nogain(i,:), fit_nogain(i).dev, fit_nogain(i).stats] = ...
            glmfit(X,y,'poisson');
        
    end

    % compute x
    preds_nogain = X * coeffs_nogain(i,2:end)'; % no intercept!
    
    % add stuff to struct
    fit_nogain(i).coeffs = coeffs_nogain(i,:);
    fit_nogain(i).pred = preds_nogain;


    %% second step: fit the gain
    
    % make dm for the test trials
    p_test = X0(ops.include & testI,:) * coeffs_nogain(i,2:end)';
    [C_test,C1_test] = dummycode(C0(ops.include & testI,:),cc0(ops.include & testI),ops);
    dm_test = [p_test, p_test.*C_test, C1_test];
    y_test = y0(ops.include & testI);
    
    % generate contrast and stimulus predictors using STRF prediction
    XC = preds_nogain .* C;
    dm_gain = [preds_nogain, XC, C1];
    
    % fit
    if ops.alpha(2) >= 0 && ops.alpha(2) <= 1
        %fprintf('Step 2: fitting gain using glmnet... ');
        
        % use regularization
        opts = glmnetSet;
        opts.alpha = ops.alpha(2);
        nfolds = 10;
        fit_gain(i) = cvglmnet(dm_gain,y,'poisson',opts,[],nfolds,[],false);
        coeffs_gain(i,:) = cvglmnetPredict(fit_gain(i),[],'lambda_1se','coefficients');
        fit_gain(i).pred_train = cvglmnetPredict(fit_gain(i),dm_gain,'lambda_1se','response');
        fit_gain(i).pred_test = cvglmnetPredict(fit_gain(i),dm_test,'lambda_1se','response');
        
    else
        %fprintf('Step 2: fitting gain using glmfit... ');
        
        % unregularized
        [coeffs_gain(i,:), fit_gain(i).dev, fit_gain(i).stats] = ...
            glmfit(dm_gain,y,'poisson');
        fit_gain(i).pred_train = glmval(coeffs_gain(i,:)',dm_gain,'log');
        fit_gain(i).pred_test = glmval(coeffs_gain(i,:)',dm_test,'log');

    end

    % correlation for training
    fit_gain(i).corr = corr(fit_gain(i).pred_train,y);
    fit_gain(i).r2 = fit_gain(i).corr.^2;
    
    % correlation for testing
    fit_gain(i).corr_test = corr(fit_gain(i).pred_test,y_test);
    fit_gain(i).r2_test = fit_gain(i).corr_test.^2;
    
    % add coeffs
    fit_gain(i).coeffs = coeffs_gain(i,:);
    
    % prune fit stats if cross-validating
    if ops.cvfolds > 1
        fit_gain(i).stats = [];
        fit_nogain(i).stats = [];
    end
    
    fprintf('%d ',i);

end
toc;


%% post processing
fprintf('Post processing data... '); tic;

% average model coefficients
beta0 = mean(coeffs_gain(:,1),1);
beta1 = mean(coeffs_gain(:,2),1);
beta2 = mean(coeffs_gain(:,3:3+size(XC,2)-1),1)';
beta3 = mean(coeffs_gain(:,3+size(XC,2):end),1)';

% old gain index (when C isn't masked)
%w_norm = (beta1+C*beta2) / (beta1+sum(beta2));

% gain index (ratio of gradient of changing contrast model to
% gradient of a model with constant contrast)
[C,C1] = dummycode(C0(ops.include,:),cc0(ops.include),ops);
C_new = C;
C_new(C ~= 0) = 1; % constant contrast model
w_norm = (beta1+C*beta2) ./ (beta1+C_new*beta2);
w = nan(size(ops.include));
w(ops.include) = w_norm;

% average strf model
res.strf = reshape(mean(coeffs_nogain(:,2:end),1),length(ops.f),[]);
res.strf_fit.coeffs = mean(coeffs_nogain,1);
res.strf_fit.pred = X0 *  res.strf_fit.coeffs(2:end)';
res.strf_fit.corr = corr(res.strf_fit.pred,y0);
res.strf_fit.r2 = res.strf_fit.corr.^2;

% average gain model
dm = [res.strf_fit.pred, res.strf_fit.pred .* C, C1];
res.gain_fit.coeffs = [beta0; beta1; beta2; beta3];
res.gain_fit.pred = glmval(res.gain_fit.coeffs,dm,'log');
res.gain_fit.corr = corr(res.gain_fit.pred,y0);
res.gain_fit.r2 = res.gain_fit.corr.^2;

% cross-validation results
res.strf_fit_cv = fit_nogain;
res.gain_fit_cv = fit_gain;

% mains
res.beta0 = beta0;
res.beta1 = beta1;
res.beta2 = beta2;
res.beta3 = beta3;
res.w = w;

% options
ops.cv_trials = cv_trials;
ops.cv_sample = cv_sample;

toc;
