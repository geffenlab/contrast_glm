clear all; close all;
addpath(genpath('~/chris-lab/code_general/'));
addpath(genpath('~/gits/kilosort-code'));
addpath(genpath('./_functions/'));


%% SETUP
ops.root1 = '~/data/kilosort';
ops.mouseList = {'K184'};
ops.mouseLabel = {'none'};
ops.datDir = '~/data/gain_opto';
ops.template_dir = 'blockTemplates.mat';
ops.drift_correct = true; % drift correct to match spikes to stim

% for each mouse, generate a path list of sorted data
ops.pathList = getSortedSessions(ops.root1,ops.mouseList);

% load and format data
if ~exist(fullfile(ops.datDir,'k184.mat'))
    [data.spikeData, data.sessionData] = catSortedCells(ops.pathList,ops);
    save(fullfile(ops.datDir,'k184.mat'),'data','ops');
else
    load(fullfile(ops.datDir,'k184.mat'));
end

% neuron loop
lastSession = 'firstSession';
for i = 1:size(data.spikeData,2)
    
    fprintf('Cell %d/%d... ',i,size(data.spikeData,2));
    
    testCell = i;
    id2use = data.spikeData(i).cellID;
    u = data.spikeData(testCell);

    %% update session data
    if ~strcmp(lastSession,u.sessionID)
        
        fprintf('\tsetting up...'); tic;
        % match to session data
        sI = contains({data.sessionData.sessID},u.sessionID);
        s = data.sessionData(sI);

        % first, format spectrogram and extract stim info
        bI = contains({s.block.name},'contrastAlt');
        spec = s.block(bI).stimInfo.DB;
        order = s.block(bI).stimInfo.allCond;
        stimInfo = s.block(bI).stimInfo.params;
        nreps = s.block(bI).nreps;
        
        lastSession = u.sessionID;
        toc;
        
    end
    u.spikes = u.spikes{bI};
    
    fn = fullfile(ops.datDir,[u.cellID '.mat']);
    save(fn,'u','s','spec','order','stimInfo','nreps');
    
end