% objectives:
% 1) given a STRF, a stimulus and a set of NL parameters, generate spikes
% 2) compute a trial averaged PSTH, and make a function for fitting
% the parameters in each contrast and the adaptation time
addpath(genpath('./_functions/'))

% stimulus generation
ntrials = 100;
trial_length = 2;
contrast = [1 3];
T = .025;
fs = 1/T;
freqs = 1:32;

S = []; C = [];
for i = 1:ntrials
    chords = [randn(length(freqs),trial_length*fs)*contrast(1) ...
              randn(length(freqs),trial_length*fs)*contrast(2)];
    S = [S ...
         chords];
    C = [C ...
         repmat(contrast(1),1,trial_length*fs) ...
         repmat(contrast(2),1,trial_length*fs)];
end


% $$$ % strf generation
% $$$ t = 1:13;
% $$$ f = freqs;
% $$$ [X1,X2] = meshgrid(t,f);
% $$$ X = [X1(:) X2(:)];
% $$$ mu = [2 15];
% $$$ Sigma = [0.7 0.1; 0.1 0.7];
% $$$ y = mvnpdf(X,mu,Sigma);
% $$$ strf = reshape(y,length(f),length(t));
% $$$ imagesc(strf)
% $$$ 
% $$$ % convolve strf with stim to get linear output
% $$$ X = makeDesignMatrix(S,length(t));
% $$$ yl = X' * strf(:);

% load test neuron
c = 15;
data = load('acute_res_spline-3-3_2tau.mat');
mdl = data.ops.ln_model;
p = 