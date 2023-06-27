clear all; close all;
addpath(genpath('./_functions/'));

filepath = '/Volumes/chris_data/timemachine_chris_imac/data/gain_opto';
c = 'CA131_2021-08-26_14-46-26_027_093_mu';
file = fullfile(filepath, [c '.mat']);

load(file);

bI = contains({s.block.name}, 'contrastALT');
spec = s.block(bI).stimInfo.DB;
order = s.block(bI).stimInfo.allCond;
stimInfo = s.block(bI).stimInfo.params;
nreps = s.block(bI).nreps;
u.spikes = u.spikes{bI};
s.block = s.block(bI);

fn = fullfile('~/chris-lab/contrast_glm/_data', [c '.mat']);

save(fn,'u','s','spec','order','stimInfo','nreps');

