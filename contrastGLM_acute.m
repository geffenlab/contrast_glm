function contrastGLM_acute(infile,outfile,w,gain_w,alpha,dummycode,cvfolds,spline_basis)
%
%% function contrastGLM_acute(infile,outfile,w,gain_w,alpha,dummycode,cvfolds,spline_basis)
%
% This function fits a GLM to the neural data in [infile], and outputs the
% fit results to the file in [outfile]. This form of the GLM is improved from
% version 2.0. Improvements are:
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
%
% inputs:
%  infile/outfile: paths to the data input and results output
%  w: window for the strf (in seconds, default .3)
%  gain_w: window for gain estimates (in seconds, default 1)
% 
% NOTES:
%  - this glm cannot be regularized... this results in bias in predicting gain control

% to compile
if false
    %%
    root = '/cbica/home/angelonc/comp_space/contrast_glm';
    addpath(genpath(fullfile(root,'_functions/')));
    disp('Compiling contrastGLM_acute.m...'); tic;
    mcc -m contrastGLM_acute.m -R -singleCompThread -R -nodisplay -R -nojvm
    cmd = sprintf('sed -i ''s|exe_dir=`dirname "$0"`|exe_dir=%s|'' ./run_contrastGLM_acute.sh',root);
    system(cmd); toc;
end



%% setup
% get host information
[~,ops.hostname] = system('hostname');

% set the seed randomly
ops.seed = rng('shuffle');

% start timer
startTime = datestr(now);
t0 = tic; fprintf('<<< Running contrastGLM_acute v.3.0 -- %s >>>\nhost: %s\n',...
                  datestr(now),ops.hostname);

% if on CUBIC login or lab computer, addpath
if contains(ops.hostname,'login') || contains(ops.hostname,'gefen')
    addpath(genpath('./_functions/'));
end

% set default file
default_input = '/cbica/home/angelonc/comp_space/_data/K184_2021-04-21_13-48-12_009_025_mu.mat';
if contains(ops.hostname,'gefen')
    default_input = '~/data/gain_opto/K184_2021-04-21_13-48-12_009_025_mu.mat';
end



%% argument handling
if nargin < 1
    infile = default_input;
end
if nargin < 2
    outfile = 'test.mat';
end
if nargin < 3
    w = .3;
end
if nargin < 4
    gain_w = 1;
end
if nargin < 5
    alpha = -1;
end
if nargin < 6
    dummycode = 1;
end
if nargin < 7
    cvfolds = 10;
end
if nargin < 8
    spline_basis = [];
end
ops.infile = infile;
ops.outfile = outfile;
ops.w = w;
ops.gain_w = gain_w;
ops.alpha = alpha;
ops.dummycode = dummycode;
ops.cvfolds = cvfolds;
ops.spline_basis = spline_basis;

% check numeric inputs for string values
inputs2check = {'w','gain_w','alpha','dummycode','cvfolds','spline_basis'};
for i = 1:length(inputs2check)
    if isstr(ops.(inputs2check{i}))
        ops.(inputs2check{i}) = str2num(ops.(inputs2check{i}));
    end
end

% duplicate alpha argument if only one is provided
if numel(ops.alpha) < 2
    ops.alpha = repmat(ops.alpha,2,1);
end



%% load the data for this neuron
fprintf('Loading data... ');
load(ops.infile); t1 = toc(t0); toc(t0);



%% parameters
ops.period = stimInfo.chordDuration;
ops.fs = 1/ops.period;
ops.t = 0:ops.period:ops.w-ops.period;
ops.gain_lags = ops.gain_w * ops.fs;
ops.f = stimInfo.freqs;
ops.stimInfo = stimInfo;
ops.contrast = stimInfo.sd;
ops.mean_contrast = 1/((ops.contrast(1)+ops.contrast(2))/(2*ops.contrast(1)*ops.contrast(2)));
disp(ops);

% add block length/chord length
ops.blockLength = stimInfo.baseNoiseD;

% shifting parameter for contrast mask
ops.shift = ops.stimInfo.baseNoiseD/2 * ops.fs;


%% stimulus formatting
fprintf('Setting up the GLM... '); tic;

% duplicate spectrogram and order
spec_r = repmat(spec',nreps,1) - mean(spec(:));
order_r = repmat(order,nreps,1);

% index the contrast,laser,scene for the whole timecourse
contrast = cleanResample(order_r(:,1)',stimInfo.baseNoiseD,stimInfo.chordDuration)';
scene = cleanResample(order_r(:,2)',stimInfo.baseNoiseD,stimInfo.chordDuration)';

% set up contrast with the true values, then normalized by its
% harmonic mean
c = contrast;
cc = ops.mean_contrast ./ c;

% make a stimulus design matrix and contrast design matrix
X = lagDesignMatrix(spec_r,length(ops.t))';

% include all times and add trial order
ops.include = true(size(c));
ops.order_r = order_r;

% get spikes and make a histogram
spikes = u.spikes;
edges = 0:stimInfo.chordDuration:length(spec_r) * stimInfo.chordDuration;
y = histcounts(spikes,edges)';
% time it
t(1) = toc(t0); toc;



%% fit two step model
if isempty(ops.spline_basis)
    [res,ops] = fitGLM(X,y,cc,ops);
else
    [res,ops] = fitGLM_bspline(X,y,cc,ops);
end
t(2) = toc(t0);


%% fit linear-nonlinear models
[res,ops] = fitLN(res,y,cc,ops);
t(3) = toc(t0);



%% compute noise ratio for spikes
y2 = reshape(y,3*ops.fs,[])';
[NR,uPattI] = responsePower(y2(:,ops.fs:end),order_r);
nr(:,1) = NR(uPattI(:,1) == ops.contrast(1));
nr(:,2) = NR(uPattI(:,1) == ops.contrast(2));

% results structure
res.y = y;
res.c = c;
res.scene = scene;
res.nr = nr;
res.cellID = u.cellID;
res.sessionID = u.sessionID;

% stop time
stopTime = datestr(now);

% add options and run time to ops struct
ops.runTime = t;
ops.timeStamps = {startTime stopTime};
ops.dataFile = ops.infile;

% save out the data
save(ops.outfile,'res','ops','u','s','spec','order','nreps');
toc(t0);

% plot_acute_neuron('test.mat')

