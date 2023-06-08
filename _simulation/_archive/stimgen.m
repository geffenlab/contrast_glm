function [stim,c] = stimgen(ops)

%% stim = stimgen(ops)
%
% makes stimuli of alternating contrasts from uniform distributions
%
% ops.mu: vector of stimulus mean
% ops.sd: vector of contrast values
% ops.blockLength: time (s) per contrast block
% ops.f: number of frequency bins
% ops.exemplars: number of frozen noise patterns
% ops.ntrials: number of repeats of each exemplar
% ops.pad: time (s) to pad each trial
% ops.fs: sample rate

% preallocate
tsamps = ops.blockLength*ops.fs*length(ops.sd)+ops.pad*ops.fs; % samples per trial
stim = zeros(tsamps * ops.exemplars * ops.ntrials,length(ops.f));
c = zeros(size(stim,1),1);

% generate stimulus
for i = 1:ops.exemplars
    
    db = [];
    for j = 1:length(ops.mu)
        db = [db unifrnd(ops.mu(j)-ops.sd(j),...
                         ops.mu(j)+ops.sd(j),...
                         length(ops.f),...
                         ops.blockLength*ops.fs)];
    end
    
    db = [db, zeros(length(ops.f),ops.pad*ops.fs)];
    s = size(db,2);
    I = (i-1)*s*ops.ntrials+1:i*s*ops.ntrials;
    stim(I,:) = repmat(db',ops.ntrials,1);
    c(I) = repmat([ones(1,ops.blockLength*ops.fs)*ops.sd(1) ...
                   ones(1,ops.blockLength*ops.fs)*ops.sd(2) ...
                   zeros(1,ops.pad*ops.fs)],...
                  1,ops.ntrials);

end

