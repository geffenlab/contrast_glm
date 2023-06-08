function [stim,St,y,contrast,contrast_inst,trialI,stimLength,ops] = buildStimAndSpikesv2(d,ops,spec,stimInfo,lvl)

%% function [stim,St,y,contrast,contrast_inst,trialI,stimLength,ops] = buildStimAndSpikesv2(d,ops,spec,stimInfo)
%
% This formats trial-wise and neural data for GLM fitting.
%
% INPUTS:
%  d = neural data structure
%  ops = fitting options
%  spec = stimulus spectrogram cell array
%  stimInfo = parameters for the stimulus set
%
% OUTPUTS:
%  stim = centered stimulus
%  St = original stimulus
%  y = spikes
%  contrast = ideal contrast, relative to the mean (ideal == width of
%             uniform dist. for stimulus creation)
%  contrast_inst = instantaneous contrast, relative to mean
%  c = contrast index (0 for no stim, -1 high contrast, 1 low contrast)
%  t = target index
%  stimLength = number of samples per trial
%  ops = new ops struct output

if nargin < 5
    lvl = false;
end

ttype = d.trialType;
for i = 1:length(ttype)
    
    if ttype(i,1) > 0
        ind = ttype(i,1);
    else
        if lvl
            ind = d.level(i);
        else
            ind = 1;
        end
    end
    
    stimIndex(i,:) = [(ttype(i,1)>0)+1,ttype(i,2),ttype(i,3),ind];
    
    % stim length per trial
    stimLength(i) = length(spec{stimIndex(i,1),stimIndex(i,2),stimIndex(i,3),stimIndex(i,4)});
    
end

% preallocate stimulus traces
noiseTrials = ttype(:,1) < 1;
tlen = length(stimLength)*ops.pad*ops.fs + sum(stimLength);
S = zeros(tlen,length(ops.f));
c = zeros(tlen,1);
t = zeros(tlen,1);
y = zeros(tlen,1);
trialI = zeros(tlen,8);
St = S;

% contrast order
if diff(stimInfo.sd) > 0
    % low to high
    cord = [1 -1];
    col = {'b','r'};
else
    % high to low
    cord = [-1 1];
    col = {'r','b'};
end

% generate stimulus matrix and indices
pad = ops.pad*ops.fs;
start = 1;
lastLength = 0;
for i = 1:length(ttype)
    
    % stimulus (always just noise)
    S(start:(start+stimLength(i)-1),:) = spec{stimIndex(i,1),...
                        stimIndex(i,2),...
                        stimIndex(i,3),...
                        1}';
    
    % stimulus with targets
    St(start:(start+stimLength(i)-1),:) = spec{stimIndex(i,1),...
                        stimIndex(i,2),...
                        stimIndex(i,3),...
                        stimIndex(i,4)}';
    
    % contrast
    c(start:(start+3*ops.fs-1)) = cord(1);
    c((start+3*ops.fs):(start+stimLength(i)-1)) = cord(2);
    
    % target time
    targI = ((d.offsets(stimIndex(i,2)) + 3) * ops.fs);
    t(start+targI-1) = 1;
    tInd(i) = start+targI;
    
    %% other trial indices
    % trial ID
    trialI(start:(start+stimLength(i)+pad-1),1) = i;
    trialI(start:(start+stimLength(i)+pad-1),2) = ...
        ttype(i,1);
    trialI(start:(start+stimLength(i)+pad-1),3) = ...
        ttype(i,2);
    trialI(start:(start+stimLength(i)+pad-1),4) = ...
        ttype(i,3);
    trialI(start:(start+stimLength(i)+pad-1),5) = ...
        d.response(i);
        
    % index same periods of time in each contrast for this trial
    % for STRF fitting (1s in C1, 1s in C2)
    trialI(start:(start+stimLength(i)+pad-1),6) = 0;
    trialI(start+(2-ops.w)*ops.fs:start+(3-ops.w)*ops.fs,6) = 1;
    trialI(start+(3+ops.w)*ops.fs:start+(3+ops.w+1)*ops.fs,6) = 1;
    
    % index the entire contrast transition up to the length of this
    % trial
    ss = stimLength(i) - ops.fs*3; % n samples after switch
    trialI(start:(start+stimLength(i)+pad-1),7) = 0;
    trialI(start+(3*ops.fs-ss):start + (3*ops.fs+ss)-1,7) = 1;
    
    % index the entire stimulus minus trial onset
    trialI(start:(start+stimLength(i)+pad-1),8) = 0;
    trialI(start+(ops.pad*ops.fs):(start+stimLength(i)+length(ops.t)-1),8) = 1;
    
    % index the target sample
    trialI(start:(start+stimLength(i)+pad-1),9) = 0;
    trialI(tInd(i),9) = 1;
    
    % index just the stimulus (not the pad)
    trialI(start:(start+stimLength(i)+pad-1),10) = 0;
    trialI(start:(start+stimLength(i)-1),10) = 1;    

    
    %% spikes
    ev = d.trialEvents(i);
    lens = stimLength(i)*ops.period;
    edges = 0:ops.period:lens;
    y(start:(start+stimLength(i)-1)) = ...
        round(makePSTH(d.spikes,ev,edges)/ops.fs);
    
    % next trial
    lastLength = stimLength(i);
    start = start + lastLength + ops.pad*ops.fs;
    
    % plot
    if false
        figure(999); clf; hold on;
        plot(trialI(1:start,:));
        plot(mean(St(1:start,:),2)-50+1,'k');
        legend('location','sw');
        keyboard
    end
    
    
end

% centered stimulus
stim = St;
stim(c==0,:) = stimInfo.mu;
stim = stim - stimInfo.mu;

% compute ideal relative contrast (ie. true contrast of distribution)
ideal_contrast    = c;
contrast(c == -1) = max(ops.contrast);
contrast(c == 1)  = min(ops.contrast);
contrast(c == 0)  = ops.mean_contrast;

% compute instantaneous relative contrast (ie. contrast at each bin)
contrast_inst          = std(stim,0,2);
ops.inst_contrast(1)   = mean(contrast_inst(c==-1));
ops.inst_contrast(2)   = mean(contrast_inst(c==1));
if ops.contrast(1) > ops.contrast(2)
    ops.inst_contrast  = [max(ops.inst_contrast) min(ops.inst_contrast)];
else
    ops.inst_contrast  = [min(ops.inst_contrast) max(ops.inst_contrast)];
end
ops.mean_inst_contrast = 1/((ops.inst_contrast(1)+ops.inst_contrast(2))/(2*ops.inst_contrast(1)*ops.inst_contrast(2)));
contrast_inst(c == 0)  = ops.mean_inst_contrast;