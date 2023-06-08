function [res,ops] = fitGLM_bspline(X0,y0,cc0,ops)
    
%% function [res,ops] = fitGLM_bspline(X0,y0,cc0,ops)
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
% ops.spline_basis - [knots degree] for bspline basis set

tic;  
    
%% SETUP
% make contrast design with spline convolution
[C0,basis] = spline_convolution(cc0,ops.gain_lags, ...
                                ops.spline_basis(1), ...
                                ops.spline_basis(2));

% dummy coded contrast and original
[C,C1] = dummycode(C0,cc0,ops);

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
    

    %% first step: fit the STRF
    [coeffs_nogain(i,:), fit_nogain(i).dev, fit_nogain(i).stats] = ...
        glmfit(X0(ops.include & trainI,:),y0(ops.include & trainI),'poisson');

    % compute full stimulus prediction
    preds_nogain = X0 * coeffs_nogain(i,2:end)'; % no intercept!
    
    % add stuff to struct
    fit_nogain(i).coeffs = coeffs_nogain(i,:);
    fit_nogain(i).pred = preds_nogain(ops.include & trainI);

    %% second step: fit the gain
    
    % full, z-scored design matrix
    dm_gain = [preds_nogain preds_nogain .* C C1];
    dm_nz = nan(size(dm_gain)); dm_nz(dm_gain ~= 0) = 1;
    dm_location = mean(dm_gain.*dm_nz,1,'omitnan');
    dm_scale = std(dm_gain.*dm_nz,1,'omitnan');
    dm_gain = (dm_gain.*dm_nz-dm_location)./dm_scale;
    dm_gain(isnan(dm_gain)) = 0;
    
    % fit
    [coeffs_gain(i,:), fit_gain(i).dev, fit_gain(i).stats] = ...
        glmfit(dm_gain(ops.include & trainI,:),y0(ops.include & trainI),'poisson');
    
    % predictions
    fit_gain(i).pred_train = glmval(coeffs_gain(i,:)',dm_gain(ops.include & trainI,:),'log');
    fit_gain(i).pred_test = glmval(coeffs_gain(i,:)',dm_gain(ops.include & testI,:),'log');

    % correlation for training
    fit_gain(i).corr_train = corr(fit_gain(i).pred_train,y0(ops.include & trainI));
    fit_gain(i).r2_train = fit_gain(i).corr_train.^2;
    
    % correlation for testing
    fit_gain(i).corr_test = corr(fit_gain(i).pred_test,y0(ops.include & testI));
    fit_gain(i).r2_test = fit_gain(i).corr_test.^2;
    
    % add coeffs
    fit_gain(i).coeffs = coeffs_gain(i,:);
    
    % scale the coefficients
    b0(i,1) = coeffs_gain(i,1) - dm_location * coeffs_gain(i,2:end)';
    b1(i,1) = coeffs_gain(i,2) ./ dm_scale(1);
    b2(i,:) = coeffs_gain(i,3:3+size(C,2)-1) ./ dm_scale(2:2+size(C,2)-1);
    b3(i,:) = coeffs_gain(i,3+size(C,2):end) ./ dm_scale(2+size(C,2):end);
    
    
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
beta0 = mean(b0,1);
beta1 = mean(b1,1);
beta2 = mean(b2,1)';
beta3 = mean(b3,1)';

% old gain index (when C isn't masked)
%w_norm = (beta1+C*beta2) / (beta1+sum(beta2));

% gain index (ratio of gradient of changing contrast model to
% gradient of a model with constant contrast)
%C_new = C;
%C_new(C ~= 0) = 1; % constant contrast model
%w_norm = (beta1+C*beta2) ./ (beta1+C_new*beta2);
%w = nan(size(ops.include));
%w(ops.include) = w_norm;

% duplicate model without gain control
cs = spline_convolution(ones(size(cc0)),ops.gain_lags, ...
                        ops.spline_basis(1), ...
                        ops.spline_basis(2));
Cs = dummycode(cs,cc0,ops);
w = (beta1+C*beta2) ./ (beta1+Cs*beta2);



% average strf model
res.strf = reshape(mean(coeffs_nogain(:,2:end),1),length(ops.f),[]);
res.strf_fit.coeffs = mean(coeffs_nogain,1);
res.strf_fit.pred = X0 *  res.strf_fit.coeffs(2:end)';
res.strf_fit.corr = corr(res.strf_fit.pred,y0);
res.strf_fit.r2 = res.strf_fit.corr.^2;

% full, scaled design matrix
dm = [res.strf_fit.pred res.strf_fit.pred .* C C1];
dm_nz = nan(size(dm)); dm_nz(dm ~= 0) = 1;
dm_location = mean(dm.*dm_nz,1,'omitnan');
dm_scale = std(dm.*dm_nz,1,'omitnan');
dm = (dm.*dm_nz-dm_location)./dm_scale;
dm(isnan(dm)) = 0;

% average gain model
res.gain_fit.coeffs = [beta0; beta1; beta2; beta3];
res.gain_fit.coeffs_unscaled = mean(cat(1,fit_gain.coeffs),1)';
res.gain_fit.pred = glmval(res.gain_fit.coeffs_unscaled,dm,'log');
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

% basis set
res.basis = basis;

% options
ops.cv_trials = cv_trials;
ops.cv_sample = cv_sample;

toc;
