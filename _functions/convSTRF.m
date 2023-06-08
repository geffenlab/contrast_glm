% 1D convolution for multiple bands
function convSS = convSTRF(S,STA,mode,g,offset)

% function convSS = convSTRF(S,STA,mode,g,offset)
% computes the linear combination of STA (neural filter) and the
% stimulus (S)
%
% g is an optional gain timecourse applied to the stimulus
%
% offset will shift the result forward by 1 sample, padding with 0

if nargin < 3
    mode = 'full';
end
if nargin < 4
    g = ones(1,size(S,2));
end
if nargin < 5
    offset = false;
end

convSS = zeros(1,size(S,2));
% for each frequency
for i = 1:size(S,1)
    % compute convolution of corresponding STA band and stim band
    x = conv(S(i,:).*g, fliplr(STA(i,:)), mode);
    %x = x(floor(length(STA)/2)+1:end-floor(length(STA)/2));
    if strcmp(mode,'full')
        x = x(1:size(S,2));
        % account for conv shift by 1 sample
        if offset
            x = [0 x(1:end-1)];
        end
    end
    convSS = convSS + x;
end
