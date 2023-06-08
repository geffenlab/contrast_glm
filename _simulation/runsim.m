clear all; close all;
chrisComp = false;
[~,name] = system('hostname');
if contains(name,'oto-gefen-007.MED.UPENN.EDU')
    chrisComp = true;
    addpath(genpath('~/chris-lab/code_general/'));
end
%addpath(genpath('/Users/chris/chris-lab/code_general/fitting/glmnet_matlab/'));



%%%%%%%%%%
%% options
ops.w = 0.5;
ops.period = 0.05;%.025;
ops.fs = 1/ops.period;
ops.f = 4000*(2 .^ ([0:16]/8));
ops.t = 0:1/ops.fs:(ops.w-ops.period);
ops.mu = [0 0];
ops.sd = [1 3];
ops.exemplars = 20;%20;
ops.blockLength = 2;
ops.ntrials = 5;%20;
ops.pad = ops.w;

trial_bins = (2*ops.blockLength+ops.pad)*ops.fs;

%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the stimulus
[stim,c] = stimgen(ops);
contrast = c; %std(stim,1,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% forward model simulation
obj.base_rate = 0.1;
obj.contrast = contrast;
obj.mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));
obj.gain_control = .5;
obj.operating_point = ops.mu(1);
obj.beta = strf_gen(ops);

% run simulation
[y,l,h,g,x0] = simulate(obj,stim,contrast);
mx = mean(stim,2);


%%%%%%%%%%%%%%
%% fit the glm
% make the stim predictors
x = lagDesignMatrix(stim-obj.operating_point,length(ops.t))';

% make the contrast predictor
cc = obj.mean_contrast./contrast;
%cc(isinf(cc)) = 1;

% make the interaction predictor
xc = lagDesignMatrix((stim-obj.operating_point).*cc,length(ops.t))';

% combine
dm = [x,xc,cc];
dm_nogain = [x];

% parallel pool for fitting
%delete(gcp('nocreate'));
%p = parpool('local',3);

% % exclude first trial
% y = y(trial_bins+1:end);
% dm = dm(trial_bins+1:end,:);
% c = c(trial_bins+1:end);

% compute fit weights based on variance of neural responses
spikevar = var(reshape(y, trial_bins, []), [], 2);
weights = 1./spikevar;
weights = repmat(weights, ops.ntrials*ops.exemplars, 1);

% first step  - fit STRF
fprintf('glmnet fit... '); tic;
opts = glmnetSet;
opts.alpha = 0.95; % 0.95% elastic net, close to lasso
nfolds = 10;
fit_nogain = cvglmnet(dm_nogain,y,'poisson',opts,[],nfolds,[],false); toc;
%coeffs = glmfit(dm,y,'poisson');%,'weights',weights);
coeffs_nogain = cvglmnetPredict(fit_nogain,[],'lambda_min','coefficients');
preds_nogain = dm_nogain * coeffs_nogain(2:end);%note no intercept%cvglmnetPredict(fit_nogain,dm_nogain,'lambda_1se','link');

preds_nogain_scale = std(preds_nogain);
preds_nogain_rescaled = preds_nogain/preds_nogain_scale;

% second step - fit gain
XC = lagDesignMatrix(preds_nogain.*cc, 5)';%length(ops.t))';
C = lagDesignMatrix(cc, 5)';%length(ops.t))';
dm_gain = [preds_nogain, XC, C];
% using glmnet
%fit_gain = cvglmnet(dm_gain, y, 'poisson', opts, [], nfolds, [], false);
%coeffs_gain = cvglmnetPredict(fit_gain,[],'lambda_min','coefficients');
% using glmfit
coeffs_gain = glmfit(dm_gain, y, 'poisson');
beta1 = coeffs_gain(2);
beta2 = coeffs_gain(3:3+size(XC,2)-1);
beta3 = coeffs_gain(3+size(XC,2):end);

% z-score no-gain preds
% dm_gain = [preds_nogain_rescaled, preds_nogain_rescaled.*cc, cc];
% coeffs_gain = glmfit(dm_gain, y, 'poisson');
% beta1 = coeffs_gain(2)/preds_nogain_scale;
% beta2 = coeffs_gain(3)/preds_nogain_scale;

% using glmnet
% opts.standardize = false;
% dm_gain = [preds_nogain_rescaled, preds_nogain_rescaled.*cc, cc];
% fit_gain = cvglmnet(dm_gain, y, 'poisson', opts, [], nfolds, [], false);
% coeffs_gain = cvglmnetPredict(fit_gain,[],'lambda_min','coefficients');
% beta1 = coeffs_gain(2)/preds_nogain_scale;
% beta2 = coeffs_gain(3)/preds_nogain_scale;



%%%%%%%%%%%%%%%%%%
%% post processing

% extract fit coefficients and prediction
%coeffs = cvglmnetPredict(fit,[],'lambda_1se','coefficients');
%pred = cvglmnetPredict(fit,dm,'lambda_1se','response');
% beta0 = coeffs(1);
% beta1 = coeffs(2:2+length(ops.f)*length(ops.t)-1);
% beta2 = coeffs(2+length(ops.f)*length(ops.t):length(coeffs)-1);
% beta3 = coeffs(end);

% debug: fake beta1
%fakebeta1 = normrnd(mean(beta1), std(beta1), size(beta1));
%fakebeta1 = zeros(size(fakebeta1));

% reshape for visualization
%b1 = reshape(beta1,length(ops.f),[]);
%b2 = reshape(beta2,length(ops.f),[]);
b1 = beta1';
b2 = beta2';

% contrast design matrix
%C = lagDesignMatrix(repmat(cc,1,length(ops.f)),length(ops.t));

% gradient-based formulation
% gamma_s = beta2 .* C; % TF x n_timepoints
% delta = gamma_s - beta2;
% v = beta1+beta2;
% w_norm = 1 + sum((v).*delta, 1)/(v'*(beta1+beta2));
% ww_norm = reshape(w_norm',(ops.blockLength*2 + ops.pad)*ops.fs,[])';

% two step, scalar
%w_norm = 1 + beta2/(beta1+beta2) * (cc-1);
% two step, matrix
w_norm = (beta1+C*beta2)/(beta1+sum(beta2));
ww_norm = reshape(w_norm',(ops.blockLength*2 + ops.pad)*ops.fs,[])';

% true gain
true_gain = zeros(size(cc));
for t=1:length(cc)
    true_gamma = zeros(size(obj.beta));
    for lag=1:size(true_gamma,2)
        lagged_t = t-lag+1;
        if lagged_t<=0
            lagged_t = lagged_t + numel(cc);
        end
        true_gamma(:,lag) = obj.beta(:,lag) * cc(lagged_t);
    end
    true_gain(t) = 1 + obj.gain_control * reshape(obj.beta, [], 1)' * reshape((true_gamma(:,:)-obj.beta), [], 1) / norm(obj.beta, 'fro')^2;
end

true_gain_instantaneous = 1 + obj.gain_control * (cc-1);

true_gain = true_gain_instantaneous;

%%%%%%%%%%%
%% plotting
f1 = figure(); clf
subplot(2,4,1); hold on;
imagesc(ops.t,ops.f/1000,obj.beta); colorbar;
title('STRF (beta)'); axis tight;

subplot(2,4,2); hold on;
scatter(h(c==ops.sd(2)),y(c==ops.sd(2)),'r.');
scatter(h(c==ops.sd(1)),y(c==ops.sd(1)),'b.');
xlabel('Linear Drive');
ylabel('Firing Rate');
title('Nonlinearity')

subplot(2,4,3); hold on;
imagesc(reshape(h,trial_bins,[])'); colorbar;
plot(repmat(ops.blockLength*ops.fs,2,1),[0 ops.ntrials*ops.exemplars],'r','linewidth',2)
xlabel('Time (step)');
ylabel('Trial');
title(sprintf('Linear Drive\nc_{low}=%02.1f,c_{hi}=%02.1f',ops.sd(1),ops.sd(2)));
axis tight;

subplot(2,4,4); hold on;
imagesc(reshape(y,trial_bins,[])'); colorbar;
plot(repmat(ops.blockLength*ops.fs,2,1),[0 ops.ntrials*ops.exemplars],'r','linewidth',2)
plot(mean(reshape(y,trial_bins,[])'),'w','linewidth',1)
xlabel('Time (step)');
ylabel('Trial');
title('Poisson Counts');
axis tight;

subplot(2,4,5); hold on;
imagesc(ops.t,ops.f/1000,b1); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('{\beta}_1'); axis tight;

subplot(2,4,6); hold on;
imagesc(ops.t,ops.f/1000,b2); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('{\beta}_2'); axis tight;

subplot(2,4,7); hold on;
plot(mean(reshape(true_gain, trial_bins, [])'),'r:', 'LineWidth', 1.5);
plot(mean(ww_norm),'k', 'LineWidth', 1.5);
legend('Gain (true)', 'Gain (estimated)')
%plot(mean(ww2),'r');
title('Gain Modulation')
% lh = legend('$w = 1 + {\beta}_2 / ({\beta}_1 + {\beta}_2) * (C-1)$',...
%             'Interpreter','latex');
ylabel('$w$','Interpreter','Latex');
xlabel('Time (step)');
%lh.Position(1) = lh.Position(1) + .2;

if chrisComp
    saveFigPDF(f1,[1000 400],...
               sprintf('sim_res_fr%d_c%d-%de-1.pdf',...
                       obj.base_rate, ...
                       ops.sd(1)*10,ops.sd(2)*10));
end























%% functions


% gain function
function g = gain(obj,contrast)
    g = obj.gain_control * obj.mean_contrast./contrast +...
        (1-obj.gain_control) * 1;
end



% gain scaling function
function [lambda_scaled, this_gain] = gain_scaling(obj,lambda,contrast)
    this_gain = gain(obj,contrast);
    lambda_scaled = log(obj.base_rate) + this_gain .* (lambda-log(obj.base_rate));
end



% firing rate function
function [l,linear_drive,g,x0] = lambda(obj,x,contrast)
    % x = design matrix
    % beta = strf weights    
    x0 = convSTRF(x'-obj.operating_point,fliplr(obj.beta))';
    linear_drive = log(obj.base_rate) + x0;
    % gain works by rescaling beta. At any given lag, the beta at that lag
    % for all frequencies is rescaled depending on the contrast at that lag
%     g = gain(obj,contrast);
%     drive_with_gain = log(obj.base_rate) + convSTRF(x'-obj.operating_point,fliplr(obj.beta), 'full',g')';
    [drive_with_gain, g] = gain_scaling(obj, linear_drive, contrast);
    l = exp(drive_with_gain);
end



% simulation
function [y,l,h,g,x0] = simulate(obj,x,contrast)
    % y = spike count
    % l = spike rate
    % h = linear drive
    % g = gain
    [l,h,g,x0] = lambda(obj,x,contrast);
    y = poissrnd(l);
end



% generate a model strf
function beta = strf_gen(ops)
    % model strf as 2d gaussian
    mu = [5, 2]; % [20 2];
    sig = [0.8 0.1; 0.1 .5];
    [X1,X2] = meshgrid(1:length(ops.f),1:length(ops.t));
    X = [X1(:) X2(:)];
    yy = mvnpdf(X,mu,sig);
    realSTA = reshape(yy,length(ops.t),length(ops.f))';
    noise = .025*randn(size(realSTA,1),size(realSTA,2));
    realSTA = realSTA + noise;
    realSTA = realSTA - mean(realSTA(:));
    beta = realSTA ./ sqrt(sum(realSTA(:).^2));
end



% 1D convolution for multiple bands
function convSS = convSTRF(S,STA,mode,g)

    % function convSS = convSTRF(S,STA)
    % computes the linear combination of STA (neural filter) and the
    % stimulus (S)
    %
    % g is optional gain timecourse.
    if ~exist('mode','var')
        mode = 'full';
    end
    if nargin < 4
        g = ones(1,size(S,2));
    end
    convSS = zeros(1,size(S,2));
    % for each frequency
    for i = 1:size(S,1)
        % compute convolution of corresponding STA band and stim band
        x = conv(S(i,:).*g, fliplr(STA(i,:)), mode);
        %x = x(floor(length(STA)/2)+1:end-floor(length(STA)/2));
        if strcmp(mode,'full')
            x = x(1:size(S,2));
            % account for conv shift by 1 sample
            x = [0 x(1:end-1)];
        end
        convSS = convSS + x;
    end
end



% design matrix lag
function X = lagDesignMatrix(S,lags)

    if size(S,2) < size(S,1)
        S = S';
    end

    nfs = size(S,1);
    X = zeros(nfs*lags,size(S,2));
    for i = 1:lags
        
        rowI = (i-1)*nfs+1:i*nfs;
        X(rowI,:) = circshift(S,(i-1),2);

    end

end



function [stim,c] = stimgen(ops)

%% stim = stimgen(ops)
%
% makes stimuli of alternating contrasts from uniform distributions
%
% ops.mu: vector of stimulus mean
% ops.sd: vector of contrast values
% ops.blockLength: time (s) per contrast block
% ops.f: number of frequency bins
% ops.exemplars: number of frozen noise patterns
% ops.ntrials: number of repeats of each exemplar
% ops.pad: time (s) to pad each trial
% ops.fs: sample rate

% preallocate
    tsamps = ops.blockLength*ops.fs*length(ops.sd)+ops.pad*ops.fs; % samples per trial
    stim = zeros(tsamps * ops.exemplars * ops.ntrials,length(ops.f));
    c = zeros(size(stim,1),1);

    mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));

    % generate stimulus
    for i = 1:ops.exemplars
        
        db = [];
        for j = 1:length(ops.mu)
            db = [db unifrnd(ops.mu(j)-ops.sd(j),...
                             ops.mu(j)+ops.sd(j),...
                             length(ops.f),...
                             ops.blockLength*ops.fs)];
        end
        
        db = [db, zeros(length(ops.f),ops.pad*ops.fs)];
        s = size(db,2);
        I = (i-1)*s*ops.ntrials+1:i*s*ops.ntrials;
        stim(I,:) = repmat(db',ops.ntrials,1);
        c(I) = repmat([ones(1,ops.blockLength*ops.fs)*ops.sd(1) ...
                       ones(1,ops.blockLength*ops.fs)*ops.sd(2) ...
                       ones(1,ops.pad*ops.fs)*mean_contrast],...%EP
                      1,ops.ntrials);

    end

end




% figure code
function saveFigPDF(h,sz,savepath);
    if size(sz,1) == 2
        sz = sz';
    end

    set(groot,{'DefaultAxesXColor','DefaultAxesYColor', ...
               'DefaultAxesZColor'},{'k','k','k'})
    set(h,'PaperPositionMode','auto');         
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','points');
    if ~isempty(sz)
        set(h,'PaperSize',sz);
        set(h,'Position',[0 0 sz]);
    end
    if exist('savepath','var') & ~isempty(savepath)
        print(h,savepath,'-dpdf','-r300','-painters');
    end
end

    
