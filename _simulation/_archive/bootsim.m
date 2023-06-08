clear all; close all;
addpath(genpath('~/chris-lab/code_general/'));


%% stimulation ops
ops.w = .5;
ops.period = .025;
ops.fs = 1/ops.period;
ops.f = 4000*(2 .^ ([0:32]/8));
ops.t = 0:1/ops.fs:(ops.w-ops.period);
ops.mu = [0 0];
ops.sd = [1 3];
ops.exemplars = 20;
ops.blockLength = 2;
ops.ntrials = 20;
ops.pad = ops.w;
ops.alpha = .95;


%% neuron model
obj.base_rate = 10;
obj.mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));
obj.gain_control = 1;
obj.operating_point = ops.mu(1);
obj.contrast = [];   % this is determined in stimgen step
obj.beta = [];       % this is determined in stimgen step


%% fit
% parallel pool for fitting
parallel = false;
if parallel
    delete(gcp('nocreate'));
    p = parpool('local',3);
end
    
for i = 1:50

    fprintf('Simulation %d/%d... ',i,50); tic;
    res(i) = runsim(ops,obj,parallel); toc;
    
    fh = gcf;
    fn = sprintf('sim_g%de-1_fr%d_c%d-%de-1_run%d.pdf',...
                 obj.gain_control*10,...
                 obj.base_rate,...
                 ops.sd(1)*10,ops.sd(2)*10,...
                 i);
    saveFigPDF(fh,[1000 600],fn);
    
end

fn = sprintf('_sim_g%de-1_fr%d_c%d-%de-1.mat',...
             obj.gain_control*10,...
             obj.base_rate,...
             ops.sd(1)*10,ops.sd(2)*10);
save(fn,'res','ops','obj');

