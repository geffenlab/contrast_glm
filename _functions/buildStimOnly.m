function [S,St,cc,t,c,trialI,stimLength] = buildStimOnly(ops,spec,stimInfo)

ttype = stimInfo.trialType;
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
trialI = zeros(tlen,4);
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
    
    % other trial indices
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),1) = i;
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),2) = ...
        ttype(i,1);
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),3) = ...
        ttype(i,2);
    trialI(start+1:(start+stimLength(i)+length(ops.t)-1),4) = ...
        ttype(i,3);  
    
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