addpath(genpath('./_functions/'))

% get example cell
ex = 'K184_2021-04-21_13-48-12_052_181_mu';
filePath = './_data';
file = fullfile(filePath, ex);

% options
outfile = fullfile(filePath, 'test.mat');
w = 0.3;
gw = 1.0;
alpha = 1;
dc = 2;
cv = 10;
spline = [3,3];

% fit cell
contrastGLM_acute(file, outfile, w, gw, alpha, dc, cv, spline)

% plot results
plot_acute_neuron(outfile)