function contrastGLM2_sim(outfile,a,gc,alpha)
    
%% function contrastGLM2_sim(outfile,a,gc,alpha)
%
% This function simulates a model neuron that adjusts its gain
% according to stimulus contrast. After simulating the neuron's
% spiking characteristics, we attempt to recover the neurons gain
% using a GLM.
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
%         b3 = operating point that depends on contrast
% 

% to compile:
if false
    addpath(genpath('~/comp_space/_code/'));
    mcc -m contrastGLM2_sim.m -R -singleCompThread -R -nodisplay -R -nojvm
end

% variable handling
if ~exist('outfile','var') || isempty(outfile)
    outfile = 'test.mat';
end

if ~exist('a','var') || isempty(a)
    a = 10;
end

if ~exist('gc','var') || isempty(gc)
    gc = 1;
end

if ~exist('alpha','var') || isempty(alpha)
    alpha = 0;
end

% convert string inputs to numbers
if ischar(a)
    a = str2num(a);
end
if ischar(gc)
    gc = str2num(gc);
end
if ischar(alpha)
    alpha = str2num(alpha);
end

fprintf('Output file: %s\nBase FR: %d\nGain control: %02.1f\nAlpha: %03.2f\n',...
    outfile,a,gc,alpha);

t0 = tic;

% set the seed randomly
seed = rng('shuffle');


%%%%%%%%%%%
%% stimulus
ops.w = .5;
ops.period = .025;
ops.fs = 1/ops.period;
ops.f = 4000*(2 .^ ([0:32]/8));
ops.t = 0:1/ops.fs:(ops.w-ops.period);
ops.mu = 30;
ops.sd = [1 3];
ops.exemplars = 20;
ops.blockLength = 2;
ops.ntrials = 20;
ops.pad = ops.w;
ops.alpha = alpha;

fprintf('Generating stimulus... ');

% make the stimulus
[stim,c] = stimgen(ops);
contrast = c; %std(stim,1,2);
t1 = toc(t0); toc(t0);

%%%%%%%%%%%%%%%
%% neuron model
obj.base_rate = a;
obj.mean_contrast = ops;
obj.gain_control = gc;
obj.operating_point_offset = 0;
obj.contrast = contrast;
obj.beta = strf_gen(ops);
obj.x0 = ops.mu;
obj.operating_point = obj.x0 + obj.operating_point_offset;

fprintf('Simulating neuron... ');

% run model neuron
[y,l,h,g,x0] = simulate(obj,stim,contrast);
neuron.ops = obj;
neuron.y = y;
neuron.l = l;
neuron.h = h;
neuron.g = g;
neuron.x0 = x0;
t2 = toc(t0); toc(t0);


%%%%%%%%%%%%%%
%% fit the glm
fprintf('Fitting GLM... ');

% set stim ITI to mean stim value 
stim(stim==0) = obj.x0;

% make the stim predictors
x = lagDesignMatrix(stim-obj.operating_point,length(ops.t))';

% make the contrast predictor
cc = obj.mean_contrast./contrast;
cc(isinf(cc)) = 1;

% make the interaction predictor
xc = lagDesignMatrix((stim-obj.operating_point).*cc,length(ops.t))';

% combine
dm = [x,xc,cc];

% fit
opts = glmnetSet;
opts.alpha = ops.alpha; % 0 = ridge, 1 = lasso
fit = cvglmnet(dm,y,'poisson',opts,[],[],[]);
t3 = toc(t0); toc(t0);



%%%%%%%%%%%%%%%%%%
%% post processing
fprintf('Post processing data... ');

% extract fit coefficients and prediction
coeffs = cvglmnetPredict(fit,[],'lambda_min','coefficients');
pred = cvglmnetPredict(fit,dm,'lambda_min','response');
beta0 = coeffs(1);
beta1 = coeffs(2:2+length(ops.f)*length(ops.t)-1);
beta2 = coeffs(2+length(ops.f)*length(ops.t):length(coeffs)-1);
beta3 = coeffs(end);

% reshape for visualization
b1 = reshape(beta1,length(ops.f),[]);
b2 = reshape(beta2,length(ops.f),[]);

% contrast design matrix
C = lagDesignMatrix(repmat(cc,1,length(ops.f)),length(ops.t));

% gain index (not sure if this is correct)
% w = 1 + beta2'./(beta1' + beta2') * (C-1);

% correct averaging (matrix sum, divided by number of filter
% parameters to get the correct average value):
strf = beta1 + beta2;
strf(strf==0) = 1e-5;
gk = beta2' ./ strf';
w = 1 + (gk*(C-1)) / ...
     (length(ops.f)*length(ops.t));
wm = reshape(w',(ops.blockLength*2 + ops.pad)*ops.fs,[])';

t4 = toc(t0); toc(t0);


% results structure
res.fit = fit;
res.beta0 = beta0;
res.beta1 = beta1;
res.beta2 = beta2;
res.beta3 = beta3;
res.b1 = b1;
res.b2 = b2;
res.gk = gk;
res.wm = wm;
res.neuron = neuron;
res.t = [t1 t2 t3 t4];
res.seed = seed;

% save
save(outfile,'res','ops','obj')

fprintf('File saved... done :)\n');


    
