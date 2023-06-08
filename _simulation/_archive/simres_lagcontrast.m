clear all; close all;
addpath(genpath('~/chris-lab/code_general/'))

% run analysis for different parameters
gc = [0 .5 1];
simtype = 'sim3'
alpha = 0.95;

% load first file
f = dir(fullfile('_res',sprintf('%s_*g%02d*alpha%02d*.mat',simtype,gc(1)*10,alpha*100)));
load(fullfile(f(1).folder,f(1).name))

% preallocate
nt = length(ops.t); nf = length(ops.f);
w = zeros(size(res.wm,2),length(f),length(gc));
g = w;
err = zeros(length(f),length(gc));
b0 = zeros(length(gc),length(f));
b1 = zeros(length(gc),length(f),nt*nf);
b2 = b1;
b3 = zeros(length(gc),length(f),nt);

% for each gain control value
for i = 1:length(gc)
    
    % build file list
    files = dir(fullfile('_res',sprintf('%s_*g%02d*alpha%02d*.mat',simtype,gc(i)*10,alpha*100)));
    
    for j = 1:length(files)
        
        load(fullfile(files(j).folder,files(j).name));
        
        % compile time course, true gain, compute error, range etc
        w(:,j,i) = mean(res.wm,1);
        g(:,j,i) = res.neuron.g(1:180);
        err(j,i) = immse(squeeze(w(:,j,i)),squeeze(g(:,j,i)));
        b0(i,j) = res.beta0;
        b1(i,j,:) = res.beta1;
        b2(i,j,:) = res.beta2;
        b3(i,j,:) = res.beta3;
                
        disp([i j]);
        
    end
    
end

f1 = figure(1); clf;

% w error
subplot(4,3,1)
cc = repmat([.66 .33 0],3,1)';
hold on;
for i = 1:length(gc)
    ph(i) = plot(mean(g(:,:,i),2),'--','Color',cc(i,:),'LineWidth',1);
    patchErrorBars(1:180,squeeze(w(:,:,i))',cc(i,:),'prctile');
end
legend(ph,'g=0.0','g=0.5','g=1.0');
ylabel('w,g');
xlabel('Time (step)'); plotPrefs;
title(sprintf('alpha = %3.2f (mean)',alpha));
plotPrefs;

% beta0 (baseline)
subplot(4,3,2); hold on;
plot([.5 3.5],repmat(log(obj.base_rate),1,2),'k--');
for i = 1:length(gc)
    errorbar(i,mean(b0(i,:)),...
             mean(b0(i,:)) - prctile(b0(i,:),2.7),...
             mean(b0(i,:)) - prctile(b0(i,:),97.5),...
             'color',cc(i,:),'linewidth',1,'marker','.','markersize',20);
end
xlim([.5 3.5]);
title('Baseline Rate');
set(gca,'xtick',1:length(gc));
set(gca,'xticklabels',num2str(gc'));
xlabel('Gain Control'); ylabel('$\beta_0$','interpreter', ...
                               'latex');
plotPrefs;

% beta 3 (contrast offset)
subplot(4,3,3); hold on;
for i = 1:length(gc)
    patchErrorBars(1:20,squeeze(b3(i,:,:)),cc(i,:),'prctile');
end
xlabel('Time (steps)');
ylabel('$\beta_3$','interpreter','latex');
title('Contrast Filter');
axis tight; plotPrefs;

% similarity of w across runs at high gain
subplot(4,3,4)
D = log(pdist(w(:,:,3)','cosine'));
Z = squareform(D);
imagesc(Z); colorbar; plotPrefs;
title('Dissimilarity of w');

for i = 1:length(gc)
    
    % beta1 + beta2
    bb1 = reshape(squeeze(mean(b1(i,:,:),2)),length(ops.f),[]);
    bb2 = reshape(squeeze(mean(b2(i,:,:),2)),length(ops.f),[]);
    
    % strf/gain kernel
    strf = bb1 + bb2;
    gk = bb2 ./ (bb1 + bb2);
    
    if i == 1
        clim1 = [min(strf(:)) max(strf(:))];
        clim2 = [min(gk(:)) max(gk(:))];
    end
    subplot(4,3,6+i);
    imagesc(ops.t,ops.f/1000,strf,clim1); colorbar;
    xlabel('Time (s)'); ylabel('Freq. (kHz)');
    title(sprintf('STRF (contrast %02.1f)',gc(i)));
    axis square; plotPrefs;
    
    subplot(4,3,9+i);
    imagesc(ops.t,ops.f/1000,gk,clim2); colorbar;
    xlabel('Time (s)'); ylabel('Freq. (kHz)');
    title(sprintf('Gain Kernel (contrast %02.1f)',gc(i)));
    axis square; plotPrefs;

end




saveFigPDF(f1,[900 900],sprintf('./_plots/_%s_res_laggedcontrast_alpha%03d.pdf', ...
                                simtype,alpha*100))