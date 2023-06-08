function contrastGLM_sim(varargin)
    
%% function contrastGLM_sim(varargin)
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
%         b3 = operating point that depends on lagged contrast
% 

% to compile:
if false
    addpath(genpath('../_functions/'));
    root = '/cbica/home/angelonc/comp_space/contrast_glm/_simulation';
    disp('Compiling contrastGLM_sim.m...'); tic;
    mcc -m contrastGLM_sim.m -R -singleCompThread -R -nodisplay -R -nojvm
    cmd = sprintf('sed -i ''s|exe_dir=`dirname "$0"`|exe_dir=%s|'' ./run_contrastGLM_sim.sh',root);
    system(cmd); toc;
end

% get host information
[~,ops.hostname] = system('hostname');

% if on CUBIC login or lab computer, addpath
if contains(ops.hostname,'login') || contains(ops.hostname,'gefen')
    addpath(genpath('../_functions/'));
end


%% set defaults
% addpath(genpath('../_functions/'));
% contrastGLM_sim('test.mat',30,[1 3],100,5,.3,0,40,-1,1,10,[7 3],0.1,1,[.05 .5])
defaults = {'test.mat',30,[1 3],100,5,.3,0,40,-1,1,10,[7 3],0.1,1,[.05 .5]};
varnames = {'outfile','mu','sd','exemplars','ntrials','w','pad','gain_lags',...
            'alpha','dummycode','cvfolds','spline_basis','base_rate','gain_control','tau'};




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

% set stim options based on defaults
for i = 1:12
    ops.(varnames{i}) = defaults{i};
end

% set neuron options
for i = 13:15
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
ops.shift = ops.fs;

% neuron model options
obj.mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));
obj.operating_point_offset = 0;
obj.expfun = @(p,x)(p(2) + (p(1) - p(2)) .* exp(-p(3).*x));

% add mean contrast to ops
ops.mean_contrast = obj.mean_contrast;

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
    ops.include = ones((ops.blockLength*2+ops.pad)*ops.fs*ops.ntrials*ops.exemplars,1);
end
ops.include = logical(ops.include);

% duplicate alpha if only one specified
if numel(ops.alpha) < 2
    ops.alpha = repmat(ops.alpha,1,2);
end

% create "trial order" to evenly cross-validate across contrast
ops.order_r = repmat(ops.sd',ops.exemplars*ops.ntrials,1);

% cross-validation
ops.cvfolds = 1;
    


%%%%%%%%%%%
%% stimulus
t0 = tic;

fprintf('Generating stimulus... ');

% make the stimulus
tic; [stim,c] = stimgen(ops); toc;
t(1) = toc(t0);



%%%%%%%%%%%%%%%
%% neuron model
% stimulus dependent parameters (others are previously loaded from
% paramfile)
obj.contrast = c;
obj.beta = strf_gen(ops);
obj.x0 = ops.mu;
obj.operating_point = obj.x0 + obj.operating_point_offset;

fprintf('<<< OPTIONS >>>\n');
disp(ops); disp(obj);

fprintf('Simulating neuron... ');

% run model neuron
tic; [y,l,h,g,x0] = simulate(obj,stim,ops); toc;
neuron.ops = obj;
neuron.y = y;
neuron.l = l;
neuron.h = h;
neuron.g = g;
neuron.x0 = x0;
t(2) = toc(t0);

% plot_sim_neuron(neuron,ops)


%%%%%%%%%%%%%%
%% fit the glm

% adjust silent periods so that they are at the default stimulus value
stim(stim==0) = obj.x0;

% normalize contrast to harmonic mean
cc = ops.mean_contrast ./ c;

% make the stim predictors
X = lagDesignMatrix(stim-obj.operating_point,length(ops.t))';

if isfield(ops,'spline_basis') & ~isempty(ops.spline_basis)
    % fit GLM with b-splined contrast
    [res,ops] = fitGLM_bspline(X,y,cc,ops);
else
    % fit the GLM
    [res,ops] = fitGLM(X,y,cc,ops);
end

t(3) = toc(t0);
res.t = t;


% save
save(ops.outfile,'res','ops','obj','neuron')

fprintf('File saved... done :)\n'); toc(t0);

% plot_sim_neuron('test.mat');

    
