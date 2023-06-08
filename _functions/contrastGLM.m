function contrastGLM(infile,outfile)


%% load the data for this neuron
%infile = '/cbica/home/angelonc/comp_space/contrast_glm/_data/CA102-1908021051-06-010-su.mat';
load(infile);


%% parameters
ops.period = stimInfo.chordDuration;
ops.fs = 1/ops.period;
ops.w = .5;
ops.t = 0:ops.period:ops.w;
ops.lags = length(ops.t);
ops.f = stimInfo.freqs;
ops.pad = ops.w;


%% stimulus formatting
% build full stimulus predictors
[S,St,cc,t,y,c,trialI] = buildStimAndSpikes(d,ops,spec,stimInfo);

% stimulus design matrix
X = makeDesignMatrix(St,length(ops.t));

% contrast scaled stim and design matrix
cs = St .* repmat(cc',size(S,2),1)';
CX = makeDesignMatrix(cs,length(ops.t));

% combined DM
D = [X' CX'];




%% fit the glm

% options
options = glmnetSet;
options.alpha = 0;

% cv fit with targets
fitglm = true;
fit = []; coeffs = []; pred = [];
if fitglm
    fit = cvglmnet(D,y,'poisson',options);  
    coeffs = cvglmnetCoef(fit);
    pred = glmnetPredict(fit.glmnet_fit,D,fit.lambda_1se,'response');
end

% add options to ops struct
ops.fitOptions = options;

% save out the data
save(outfile,'fit','coeffs','y','pred','ops');

% alpha = coeffs(2:(length(coeffs)-1)/2+1);
% strf = reshape(alpha,length(ops.f),[]);
% imagesc(strf)
