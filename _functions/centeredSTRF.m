function [cstrf, cf, lag] = centeredSTRF(strf,ops)

nt = length(ops.t);
nf = length(ops.f);
sf = strf;
sf(sf < 0) = nan;

% time and frequency means
pkt = mean(sf,1,'omitnan') ./ sum(mean(sf,1,'omitnan'));
pkf = mean(sf(:,1:5),2,'omitnan');

% get max frequency and lag
[~,cfi] = max(pkf); 
cf = ops.f(cfi);
lag = ops.t(pkt==max(pkt));

% index around the cf
bins = ops.bins;

% index around the cf
csfi = [cfi-bins : cfi+bins];
csi = [1:bins*2+1];
csfi = csfi(csfi <= nf & csfi > 0);
csi = csi(csfi <= nf & csfi > 0);

% extract values from strf around cf
cstrf = nan(bins*2+1,nt);
cstrf(csi,:) = strf(csfi,:);