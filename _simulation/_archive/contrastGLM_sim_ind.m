function contrastGLM_sim_ind(varargin)
    
%% function contrastGLM_sim_ind(varargin)
%
% This function simulates a model neuron that adjusts its gain
% according to stimulus contrast. After simulating the neuron's
% spiking characteristics, we attempt to recover the neurons gain
% using a GLM.
% 
% ***NOTE: in this version, the GLM is configured to separately fit
% low and high contrast ***
%
% forward model:
%  y = exp(a + g(sigma)b(x - c))
%  where y = poisson rate
%        a = log(base firing rate)
%        g = gain function that changes with contrast (sigma)
%        b = parameter determining the stimulus response (stimulus filter)
%        c = operating point, which can be contrast dependent
%
% glm:
%  y = exp(b0 + b1(x-u) + b2*(sigma_mean/sigma)*(x-u) + b3*(sigma_mean/sigma)
%   where b0 = firing rate
%         b1+b2 = the stimulus filter
%         b2/(b1+b2) = gain per filter parameter
%         b3 = operating point that depends on lagged contrast
% 

% to compile:
if false
    addpath(genpath('../_functions/'));
    root = '/cbica/home/angelonc/comp_space/contrast_glm/_simulation';
    disp('Compiling contrastGLM_sim_ind.m...'); tic;
    mcc -m contrastGLM_sim_ind.m -R -singleCompThread -R -nodisplay -R -nojvm
    cmd = sprintf('sed -i ''s|exe_dir=`dirname "$0"`|exe_dir=%s|'' ./run_contrastGLM_sim_ind.sh',root);
    system(cmd); toc;
end




%% set defaults
defaults = {'test.mat',30,[1 3],100,5,.3,0,40,-1,false,0.1,1,[.5 .05]};
varnames = {'outfile','mu','sd','exemplars','ntrials','w','pad','gain_lags',...
            'alpha','dummycode','base_rate','gain_control','tau'};




%% argument handling
% set the seed randomly
ops.seed = rng('shuffle');

% if there are no inputs, use defaults
numvarargs = length(varargin);
if numvarargs < 1
    addpath(genpath('../_functions/'));
elseif numvarargs > length(defaults)
    error('contrastGLM_sim.m:TooManyInputs','11 or less inputs');
end

% set defaults using string inputs or numeric
for i = 1:numvarargs
    if ischar(varargin{i})
        if ischar(defaults{i})
            % for string inputs
            evalcall = sprintf('defaults{%d} = ''%s'';',i,varargin{i});
        else
            % for numeric
            evalcall = sprintf('defaults{%d} = %s;',i,varargin{i});
        end
        %sdisp(evalcall);
        eval(evalcall);
    else
        defaults{i} = varargin{i};
    end
end

% set stim options based on defaults (first 7 are stim ops)
for i = 1:10
    ops.(varnames{i}) = defaults{i};
end

% set neuron options (8-11)
for i = 11:13
    obj.(varnames{i}) = defaults{i};
end




%% fixed options
% stimulus/fitting parameters
ops.period = .025;
ops.fs = 1/ops.period;
ops.f = 4000*(2 .^ ([0:32]/8));
ops.t = 0:1/ops.fs:(ops.w-ops.period);
ops.blockLength = 2;
%ops.pad = ops.w;

% neuron model options
obj.mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));
obj.operating_point_offset = 0;
obj.expfun = @(p,x)(p(2) + (p(1) - p(2)) .* exp(-p(3).*x));

% exclusion options
ops.exclude_transition = false;
if ops.exclude_transition
    % one trial
    t = ones(1,(ops.blockLength*2+ops.pad)*ops.fs);
    mx = max([ops.gain_lags,length(ops.t)]);
    t(1:mx) = 0;                             % remove onset up to max lags
    t(ops.blockLength*2*ops.fs:end) = 0;     % remove trial offset

    % repmat by number of trials
    ops.include = repmat(t,1,ops.ntrials*ops.exemplars);
else
    ops.include = ones(1,(ops.blockLength*2+ops.pad)*ops.fs*ops.ntrials*ops.exemplars);
end
ops.include = logical(ops.include);
    



%%%%%%%%%%%
%% stimulus
t0 = tic;

fprintf('Generating stimulus... ');

% make the stimulus
[stim,c] = stimgen(ops);
contrast = c; %std(stim,1,2);
t(1) = toc(t0); toc(t0);




%%%%%%%%%%%%%%%
%% neuron model
% stimulus dependent parameters (others are previously loaded from
% paramfile)
obj.contrast = contrast;
obj.beta = strf_gen(ops);
obj.x0 = ops.mu;
obj.operating_point = obj.x0 + obj.operating_point_offset;

fprintf('<<< OPTIONS >>>\n');
disp(ops); disp(obj);

fprintf('Simulating neuron... ');

% run model neuron
[y,l,h,g,x0] = simulate(obj,stim,contrast,ops);
neuron.ops = obj;
neuron.y = y;
neuron.l = l;
neuron.h = h;
neuron.g = g;
neuron.x0 = x0;
t(2) = toc(t0); toc(t0);

% plot_sim_neuron(neuron,ops)



%%%%%%%%%%%%%%
%% fit the glm

% adjust silent periods so that they are at the default stimulus value
stim(stim==0) = obj.x0;

% make the stim predictors
X = lagDesignMatrix(stim-obj.operating_point,length(ops.t))';

% contrast predictors
cc = obj.mean_contrast ./ contrast;
C = lagDesignMatrix(cc,ops.gain_lags)';

% index the predictors and spike rate
X = X(ops.include,:);
C = C(ops.include,:);
y = y(ops.include,:);

% dummy coded contrast
if ops.dummycode
    
    uc = unique(C);
    c1 = C; c2 = C;
    c1(c1 == uc(1)) = 0;
    c2(c2 == uc(2)) = 0;
    C = [c1 c2];
    
end   

% first step: fit the STRF
if ops.alpha >= 0 && ops.alpha <= 1
    fprintf('Step 1: fitting STRF using glmnet... ');
    
    % use regularization
    opts = glmnetSet;
    opts.alpha = ops.alpha;
    nfolds = 10;
    fit_nogain = cvglmnet(X,y,'poisson',opts,[],nfolds,[],false);
    coeffs_nogain = cvglmnetPredict(fit_nogain,[],'lambda_1se','coefficients');
    
else
    fprintf('Step 1: fitting STRF using glmfit... ');
    
    % unregularized
    [coeffs_nogain, fit_nogain.dev, fit_nogain.stats] = ...
        glmfit(X,y,'poisson');
    
end
preds_nogain = X * coeffs_nogain(2:end); % no intercept!

t(3) = toc(t0); toc(t0);

% second step: fit the gain
fprintf('Step 2: fitting gain using glmfit...');

% generate contrast and stimulus predictors
XC = preds_nogain .* C;
dm_gain = [preds_nogain, XC, C];

% fit contrast 1
uc = unique(cc,'stable');
c1i = cc == uc(1);
[coeffs_gain1, fit_gain1.dev, fit_gain1.stats] = ...
    glmfit(dm_gain(c1i,:),y(c1i),'poisson');
fit_gain1.pred = glmval(coeffs_gain1,dm_gain(c1i,:),'log');
fit_gain1.corr = corr(fit_gain1.pred,y(c1i));

% fit contrast 2
c2i = cc == uc(2);
[coeffs_gain2, fit_gain2.dev, fit_gain2.stats] = ...
    glmfit(dm_gain(c2i,:),y(c2i),'poisson');  
fit_gain2.pred = glmval(coeffs_gain2,dm_gain(c2i,:),'log');
fit_gain2.corr = corr(fit_gain2.pred,y(c2i));

% fit all contrasts
[coeffs_gain, fit_gain.dev, fit_gain.stats] = ...
    glmfit(dm_gain,y,'poisson');
fit_gain.pred = glmval(coeffs_gain,dm_gain,'log');
fit_gain.corr = corr(fit_gain.pred,y);
t(4) = toc(t0); toc(t0);



%%%%%%%%%%%%%%%%%%
%% post processing
fprintf('Post processing data... ');

%% complete model
% second step coefficients
beta0 = coeffs_gain(1);
beta1 = coeffs_gain(2);
beta2 = coeffs_gain(3:3+size(XC,2)-1);
beta3 = coeffs_gain(3+size(XC,2):end);

% gain index (ratio of contrast gradient to no contrast gradient)    
w_norm = (beta1+C*beta2) / (beta1+sum(beta2));
wm = nan(size(ops.include));
wm(ops.include) = w_norm;
wm = reshape(wm,[],ops.exemplars*ops.ntrials)';
t(5) = toc(t0); toc(t0);

% results structure
res.strf = coeffs_nogain;
res.strf_pred = preds_nogain;
res.strf_fit = fit_nogain;
res.gain_fit = fit_gain;
res.beta0 = beta0;
res.beta1 = beta1;
res.beta2 = beta2;
res.beta3 = beta3;
res.wm = wm;
res.t = t;


%% independent models

res.gain_fit1 = fit_gain1;
res.gain_fit2 = fit_gain2;

% contrast 1 betas and w
res.gain_fit1.beta0 = coeffs_gain1(1);
res.gain_fit1.beta1 = coeffs_gain1(2);
res.gain_fit1.beta2 = coeffs_gain1(3:3+size(XC,2)-1);
res.gain_fit1.beta3 = coeffs_gain1(3+size(XC,2):end);
res.gain_fit1.w_norm = (res.gain_fit1.beta1+C.*c1i*res.gain_fit1.beta2) / ...
    (res.gain_fit1.beta1+sum(res.gain_fit1.beta2));

% contrast 2 betas and w
res.gain_fit2.beta0 = coeffs_gain2(1);
res.gain_fit2.beta1 = coeffs_gain2(2);
res.gain_fit2.beta2 = coeffs_gain2(3:3+size(XC,2)-1);
res.gain_fit2.beta3 = coeffs_gain2(3+size(XC,2):end);
res.gain_fit2.w_norm = (res.gain_fit2.beta1+C.*c2i*res.gain_fit2.beta2) / ...
    (res.gain_fit2.beta1+sum(res.gain_fit2.beta2));

% combined w, predictions
res.ind_model.pred = nan(size(y));
res.ind_model.pred(c1i) = res.gain_fit1.pred;
res.ind_model.pred(c2i) = res.gain_fit2.pred;
res.ind_model.w_norm = nan(size(y));
res.ind_model.w_norm(c1i) = res.gain_fit1.w_norm(c1i);
res.ind_model.w_norm(c2i) = res.gain_fit2.w_norm(c2i);

% plot_sim_neuron('test.mat');

% save
save(ops.outfile,'res','ops','obj','neuron')

fprintf('File saved... done :)\n');
    
