function contrastGLM(infile,outfile,paramfile)
%
%% function contrastGLM(infile,outfile)
%
% This function fits a GLM to the neural data in [infile], and outputs the
% fit results to the file in [outfile]. This form of the GLM is improved from
% version 1.0. Improvements are:
% - stimulus is centered (x-u)
% - contrast is relative to its harmonic mean (sigma_mean/sigma)
% - inclusion of an 'operating point', a contrast-based predictor that
%   scales the firing rate depending on the contrast
%
% glm:
%  y = exp(b0 + b1(x-u) + b2*(sigma_mean/sigma)*(x-u) + b3*(sigma_mean/sigma)
%   where b0 = baseline firing rate
%         b1+b2 = the stimulus filter
%         b2/(b1+b2) = gain per filter parameter
%         b3 = operating point that depends on contrast

% to compile
if false
    %%
    root = '/cbica/home/angelonc/comp_space/contrast_glm';
    addpath(genpath(fullfile(root,'_functions/')));
    disp('Compiling contrastGLM.m...'); tic;
    mcc -m contrastGLM.m -R -singleCompThread -R -nodisplay -R -nojvm
    cmd = sprintf('sed -i ''s|exe_dir=`dirname "$0"`|exe_dir=%s|'' ./run_contrastGLM.sh',root);
    system(cmd); toc;
end

% get host information
[~,ops.hostname] = system('hostname');
    
% start timer
startTime = datestr(now);
t0 = tic; fprintf('<<< Running contrastGLMv3 -- %s >>>\nhost: %s\n',...
    datestr(now),ops.hostname);


%% load the data for this neuron
if ~exist('infile','var')
    if contains(ops.hostname,'login')
        addpath(genpath('./_functions/'));
    end
    infile = '/cbica/home/angelonc/comp_space/_data/CA104-1909181600-08-031-su.mat';
    % CA121-2004151233-25-061-su.mat';
    % CA122-2007101419-06-015-su.mat';
    % CA122-2007101419-10-022-mu.mat';
    % CA122-2007101419-03-009-su.mat';
    % CA122-2007101419-04-011-su.mat';
    % CA104-1909181600-08-031-su.mat'; % this one will fail without reg.
    outfile = 'test.mat';
    
    % run in parallel while debugging (for speed)
    %delete(gcp('nocreate'));
    %p = parpool('local',10);
    par = false;
    
else
    par = false;
    
end
fprintf('Loading data... ');
load(infile); t1 = toc(t0); toc(t0);


%% parameters
ops.period = stimInfo.chordDuration;
ops.fs = 1/ops.period;
if ~exist('paramfile','var') || isempty(paramfile) || ~exist(paramfile,'file')
    
    % default parameters
    ops.w = .25;
    ops.alpha = -1;
    ops.contrast2use = 'ideal';
    ops.gain_lags = .25 * ops.fs;
    ops.scale = false;
    ops.initializeSTRF = false;
    ops.pad = max([ops.w ops.gain_lags/ops.fs]);;

    
    % using w_norm, we do NOT want to use 'actual'....
    % ops.contrast2use = 'actual';

else
    
    % load parameters
    prms = importdata(paramfile);
    for i = 1:length(prms)
        eval(prms{i});
    end
    
end
ops.t = 0:ops.period:ops.w-ops.period;
ops.lags = length(ops.t);
ops.f = stimInfo.freqs;
ops.stimInfo = stimInfo;
ops.contrast = stimInfo.sd;
ops.mean_contrast = 1/((ops.contrast(1)+ops.contrast(2))/(2*ops.contrast(1)*ops.contrast(2)));
disp(ops);


%% stimulus formatting
fprintf('Setting up the GLM... ');

% build full stimulus predictors
[stim,~,y,contrast,~,trialI,stimLength,ops] = buildStimAndSpikesv2(d,ops,spec,stimInfo);

% select contrast to be used
if contains(ops.contrast2use,'ideal')
    cc = ops.mean_contrast ./ contrast';
else
    cc = ops.mean_inst_contrast ./ contrast_inst;
end


%% design matrix

% stimulus design matrix
X = lagDesignMatrix(stim,length(ops.t))';

% lagged contrast predictor
C = lagDesignMatrix(cc,ops.gain_lags)';

% remove data around targets, transitions, and during the pad
ops.excludeTargets = false;
ops.excludeSilentTransitions = true;
ops.include = true(size(y));
if ops.excludeTargets
    targetIndex = sum(lagDesignMatrix(trialI(:,9) & trialI(:,2)>0,ops.gain_lags))';
    ops.include = ops.include & targetIndex == 0;
elseif ops.excludeSilentTransitions
    ops.include = ops.include & trialI(:,10) == 1 & trialI(:,8);
end

% time it
t(1) = toc(t0); toc(t0);

[coeffs,dev] = glmfit(X(ops.include,:),y(ops.include),'poisson');
strf = reshape(coeffs(2:end),length(ops.f),[]);
preds_nogain = X(ops.include,:) * coeffs(2:end);

dm = [preds_nogain, preds_nogain .* C(ops.include,:), C(ops.include,:)];
[coeffs_gain,dev,stats] = glmfit(dm,y(ops.include),'poisson');
preds_gain = dm * coeffs_gain(2:end);


%% fit two step model
fits = fit2stepGLM(y,X,C,ops);
t(2) = toc(t0);


%% post processing
fprintf('Post processing data... ');

% second step coefficients
beta0 = fits.gain_fit.coeffs(1);
beta1 = fits.gain_fit.coeffs(2);
beta2 = fits.gain_fit.coeffs(3:3+size(C,2)-1);
beta3 = fits.gain_fit.coeffs(3+size(C,2):end);

% gain index (ratio of contrast gradient to no contrast gradient)  
w_norm = nan(size(ops.include));
w_norm(ops.include) = (beta1+C(ops.include,:)*beta2) / (beta1+sum(beta2));
t(3) = toc(t0); toc(t0);

% results structure
res.scale = fits.gain_fit.scale;
res.strf = fits.strf_fit.coeffs(2:end);
res.strf_fit = fits.strf_fit;
res.gain_fit = fits.gain_fit;
res.beta0 = beta0;
res.beta1 = beta1;
res.beta2 = beta2;
res.beta3 = beta3;
res.w = w_norm;
res.y = y;
res.contrast = contrast;
res.stimLength = stimLength;
res.trialI = single(trialI);

% stop time
stopTime = datestr(now);

% add options and run time to ops struct
ops.runTime = t;
ops.timeStamps = {startTime stopTime};
ops.dataFile = infile;

% plot_real_neuron('test.mat')

% save out the data
save(outfile,'res','ops','d','spec');
