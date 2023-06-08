clear all; close all;
addpath(genpath('~/chris-lab/code_general/'))
addpath(genpath('../_functions/'));

%%
% this script analyzes glm fits to neural simulations, looking at
% effects of the degree of gain control (gc) and the amount of data
% used to fit (essentially, the number of unique noise exemplars, ex)

% run analysis for different parameters
gc = [-1 -.5 0 .5 1]
ex = [5 50 100]
nt = [100 10 5]
tag = 'noregweighted';
resDir = fullfile('_res');

% load first file
f = dir(fullfile(resDir,sprintf('sim_%s_gc%03d*ex%02d*nt%03d*.mat',tag,gc(1)*10,ex(1),nt(1))));
load(fullfile(f(1).folder,f(1).name))

% preallocate
ntf = length(ops.t); nff = length(ops.f);
ns = (ops.blockLength*2 + ops.pad)*ops.fs;
w = nan(length(f),length(gc),length(ex),ns);
g = w;
strf = nan(length(f),length(gc),length(ex),length(res.strf));
b0 = nan(length(f),length(gc),length(ex));
b1 = b0;
b2 = nan(length(f),length(gc),length(ex),length(res.beta2));
b3 = b2;
y = nan(length(f),length(gc),length(ex),ns);
ld = nan(length(ex),ns*ops.ntrials*ops.exemplars);

% for each exemplar type
for i = 1:length(ex)
    
    % for each gain control value
    for j = 1:length(gc)
        
        % build file list
        files = dir(fullfile(resDir,...
                             sprintf('sim_%s*gc%03d*ex%02d*nt%03d*.mat',tag,...
                                     gc(j)*10,ex(i),nt(i))));
        
        for k = 1:length(files)
            
            load(fullfile(files(k).folder,files(k).name));
            
            % compile time course, true gain, compute error, range etc
            g(k,j,i,:) = neuron.g(1:ns);
            b0(k,j,i) = res.beta0;
            b1(k,j,i) = res.beta1;
            b2(k,j,i,:) = res.beta2;
            b3(k,j,i,:) = res.beta3;
            strf(k,j,i,:) = res.strf;
            y(k,j,i,:) = mean(reshape(neuron.y,ns,[])',1);
            
            % gain estimate (w)
            w(k,j,i,:) = mean(res.wm,1);
            
            fprintf('%d exemplars, %03.1f gain control, run %03d\n',ex(i),gc(j),k);
            
        end
        
    end
    
    ld(i,:) = neuron.h;

end


f1 = figure(1); clf;
colors = [1 0 1; ...
          .5 0 .5; ...
          .5 .5 .5; ...
          0 .5 .5; ...
          0 1 1];
for i = 1:length(gc)
    lg{i} = sprintf('g=%3.1f',gc(i));
end
nrows = 6;

%  % model strf of one neuron
%  subplot(nrows,9,3); hold on;
%  imagesc(ops.t,ops.f/1000,neuron.ops.beta); colorbar;
%  xlabel('Time (s)'); ylabel('Freq. (kHz)');
%  title('True STRF'); plotPrefs; axis tight;

for i = 1:length(ex)
    
    
    
    % linear drive for the neuron
    subplot(nrows,6,2*(i-1)+[1]); hold on;
    imagesc(reshape(ld(i,:),ns,[])'); colorbar;
    title(sprintf('%d FROZEN NOISE PATTERNS\n\nLinear Drive',ex(i))); plotPrefs;
    plotPrefs; axis tight;
    
    % strf
    subplot(nrows,6,2*(i-1)+[2]); hold on;
    ms = squeeze(mean(squeeze(strf(:,:,i,2:end)),[1 2]));
    imagesc(ops.t,ops.f/1000,reshape(ms,nff,[])); colorbar;
    xlabel('Lag (s)'); ylabel('Freq. (kHz)');
    title('STRF Estimate'); axis tight; plotPrefs;
    
    % y
    subplot(nrows,6,2*(i-1)+[7 8]); hold on;
    for ii = 1:length(gc)
        ph(ii) = patchErrorBars(1:ns,log(squeeze(y(:,ii,i,:))),colors(ii,:));
    end
    legend(ph,lg,'location','nw');
    plotPrefs; axis tight;
    ylabel('log(y)'); ylim([-3 2]);
    title('Spike Rate');

    % w
    subplot(nrows,6,2*(i-1)+[13 14]); hold on;
    for ii = 1:length(gc)
        patchErrorBars(1:ns,squeeze(w(:,ii,i,:)),colors(ii,:),'prctile');
        ph(ii) = plot(squeeze(mean(g(:,ii,i,:),1)),'--','Color', ...
                      colors(ii,:),'LineWidth',1);
    end
    plotPrefs; axis tight;
    ylabel('w,g'); ylim([0 2]);
    xlabel('Time (step)'); plotPrefs;
    title('Gain Estimate');
    
    % beta2
    subplot(nrows,6,2*(i-1)+[19 20]); hold on;
    for ii = 1:length(gc)
        ph(ii) = patchErrorBars(1:size(b2,4),squeeze(b2(:,ii,i,:)),colors(ii,:),'prctile');
        ph(ii).FaceAlpha = .1;
    end
    xlabel('Lag (step)');
    ylabel('$\beta_2$','interpreter','latex');
    title('Gain Filter');
    axis tight; plotPrefs;
    
    % beta3
    subplot(nrows,6,2*(i-1)+[25 26]); hold on;
    for ii = 1:length(gc)
        ph(ii) = patchErrorBars(1:size(b3,4),squeeze(b3(:,ii,i,:)),colors(ii,:),'prctile');
        ph(ii).FaceAlpha = .1;
    end
    ylabel('$\beta_3$','interpreter','latex');
    title('Contrast Filter');
    axis tight; plotPrefs;
    
    % beta0
    subplot(nrows,6,2*(i-1)+[31]); hold on;
    plot([.5 5.5],repmat(log(obj.base_rate),1,2),'k--');
    for ii = 1:length(gc)
        ph(ii) = errorbar(ii,mean(b0(:,ii,i)),...
                          mean(b0(:,ii,i)) - prctile(b0(:,ii,i),2.7),...
                          mean(b0(:,ii,i)) - prctile(b0(:,ii,i),97.5),...
                          'color',colors(ii,:),'linewidth',1,'marker','.','markersize',20);
    end
    xlim([.5 5.5]);
    title(sprintf('%d FROZEN NOISE PATTERNS\n\nBaseline Spike Rate',ex(i))); plotPrefs;
    set(gca,'xtick',1:length(gc));
    set(gca,'xticklabels',num2str(gc'));
    xlabel('Gain Control'); ylabel('$\beta_0$','interpreter', ...
                                   'latex');
    
    % beta1
    subplot(nrows,6,2*(i-1)+[32]); hold on;
    plot([.5 5.5],[1 1],'k--');
    for ii = 1:length(gc)
        ph(ii) = errorbar(ii,mean(b1(:,ii,i)),...
                          mean(b1(:,ii,i)) - prctile(b1(:,ii,i),2.7),...
                          mean(b1(:,ii,i)) - prctile(b1(:,ii,i),97.5),...
                          'color',colors(ii,:),'linewidth',1,'marker','.','markersize',20);
    end
    xlim([.5 5.5]); title('Stimulus Scaling');
    set(gca,'xtick',1:length(gc));
    set(gca,'xticklabels',num2str(gc'));
    xlabel('Gain Control'); ylabel('$\beta_1$','interpreter', ...
                                   'latex');
    


    
end


f2 = figure(2); clf;

% b2 lag 0
subplot(2,2,1); hold on;
edges = linspace(min(min(min(b2(:,3,:,:)))),...
                 max(max(max(b2(:,3,:,:)))),30);
histogram(b2(:,3,1,2:end),edges)
histogram(b2(:,3,2,2:end),edges)
histogram(b2(:,3,3,2:end),edges)
legend('5 examples','50 examples','100 examples');
xlabel('\beta_2[0]');
title('\beta_2 (GC = 0, lag 1:49)');
plotPrefs;

% b2 lag 2:end
subplot(2,2,3); hold on;
histogram(b2(:,3,1,1),edges)
histogram(b2(:,3,2,1),edges)
histogram(b2(:,3,3,1),edges)
legend('5 examples','50 examples','100 examples');
xlabel('\beta_2[2:end]');
title('\beta_2 (GC = 0, lag 0)');
plotPrefs;

% beta 1
subplot(2,2,2); hold on;
edges = linspace(min(min(b1(:,3,:))),...
                 max(max(b1(:,3,:))),30);
histogram(b1(:,3,1),edges)
histogram(b1(:,3,2),edges)
histogram(b1(:,3,3),edges)
legend('5 examples','50 examples','100 examples');
xlabel('\beta_1');
title('\beta_1 (GC = 0)');
plotPrefs;

% delta_w
subplot(2,2,4); hold on;
dw = w(:,:,:,ops.blockLength * ops.fs - 1) - ...
     w(:,:,:,ops.blockLength * 2 * ops.fs - 1);
edges = linspace(min(min(dw(:,3,:))),...
                 max(max(dw(:,3,:))),30);
histogram(dw(:,3,1),edges)
histogram(dw(:,3,2),edges)
histogram(dw(:,3,3),edges)
legend('5 examples','50 examples','100 examples');
xlabel('w[low] - w[high]');
title('\Deltaw (GC = 0)');
plotPrefs;


saveFigPDF(f1,[1200 1200],sprintf('./_plots/_res_2step_%s.pdf',tag));
saveFigPDF(f2,[700 450],sprintf('./_plots/_res_2step_%s_bias.pdf',tag));
