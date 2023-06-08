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

% compute mean contrast
mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));

% generate stimulus
for i = 1:ops.exemplars
    
    db = [];
    for j = 1:length(ops.sd)
        db = [db unifrnd(ops.mu-ops.sd(j),...
                         ops.mu+ops.sd(j),...
                         length(ops.f),...
                         ops.blockLength*ops.fs)];
    end
    
    % pad
    db = [db, zeros(length(ops.f),ops.pad*ops.fs)];
    
    % repeat
    s = size(db,2);
    I = (i-1)*s*ops.ntrials+1:i*s*ops.ntrials;
    stim(I,:) = repmat(db',ops.ntrials,1);
    c(I) = repmat([ones(1,ops.blockLength*ops.fs)*ops.sd(1) ...
                   ones(1,ops.blockLength*ops.fs)*ops.sd(2) ...
                   ones(1,ops.pad*ops.fs)*mean_contrast],... % CA 3/20/21 - changed from 0s
                  1,ops.ntrials);

end
