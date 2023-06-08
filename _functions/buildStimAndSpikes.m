function [S,St,cc,t,y,c,trialI,stimLength] = buildStimAndSpikes(d,ops,spec,stimInfo)

ttype = d.trialType;
for i = 1:length(ttype)
    
    if ttype(i,1) > 0
        ind = ttype(i,1);
    else
        ind = 1;
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
start = 0;
lastLength = 0;
for i = 1:length(ttype)
    
    % stimulus (always just noise)
    S(start+1:(start+stimLength(i)),:) = spec{stimIndex(i,1),...
                        stimIndex(i,2),...
                        stimIndex(i,3),...
                        1}';
    
    % stimulus with targets
    St(start+1:(start+stimLength(i)),:) = spec{stimIndex(i,1),...
                        stimIndex(i,2),...
                        stimIndex(i,3),...
                        stimIndex(i,4)}';
    
    % contrast
    c(start+1:(start+3*ops.fs)) = cord(1);
    c((start+3*ops.fs+1):(start+stimLength(i))) = cord(2);
    
    % target time
    targI = ((d.offsets(stimIndex(i,2)) + 3) * ops.fs);
    t(start+1+targI) = 1;
    tInd(i) = start+1+targI;
    
    %% other trial indices
    % trial ID
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),1) = i;
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),2) = ...
        ttype(i,1);
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),3) = ...
        ttype(i,2);
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),4) = ...
        ttype(i,3);  
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),5) = ...
        d.response(i);
        
    % index same periods of time in each contrast for this trial
    % for STRF fitting (1s in C1, 1s in C2)
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),6) = 0;
    trialI(start+1+(2-ops.w)*ops.fs:start+(3-ops.w)*ops.fs,6) = 1;
    trialI(start+1+(3+ops.w)*ops.fs:start+(3+ops.w+1)*ops.fs,6) = 1;
    
    % index the entire contrast transition up to the length of this
    % trial
    ss = stimLength(i) - ops.fs*3; % n samples after switch
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),7) = 0;
    trialI(start+1+(3*ops.fs-ss):start + (3*ops.fs+ss),7) = 1;
    
    % index the entire stimulus minus trial onset
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),8) = 0;
    trialI(start+1+(ops.w*ops.fs):(start+stimLength(i)+length(ops.t)-1),8) = 1;
    
    
    %% spikes
    ev = d.trialEvents(i);
    lens = stimLength(i)*ops.period;
    edges = 0:ops.period:lens;
    y(start+1:(start+stimLength(i))) = ...
        round(makePSTH(d.spikes,ev,edges)/ops.fs);
    
    % next trial
    lastLength = stimLength(i);
    start = start + lastLength + ops.pad*ops.fs;
    
end

% compute contrast at each timepoint (across frequencies)
cval = std(S(c~=0,:),[],2);

% normalized for regression
cn = normalize(cval.*c(c~=0),'range',[-.5 .5]);
cc = c;
cc(c~=0) = cn;