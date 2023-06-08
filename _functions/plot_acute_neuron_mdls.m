function [cor,r2] = plot_acute_neuron_mdls(res,r,u,s,ops,twin,pI,wI)

% firing rate predictions to plot
if ~exist('pI','var')
    pI = 1:size(r.pred,1);
end

% gain predictions to plot
if ~exist('wI','var')
    wI = 1:size(r.w,1);
end

% models with w (gain estimate)
w_mdls = find(contains(ops.mdls,'symm'));

% plot params
ncols = 3;
if isfield(r,'b2_kernel')
    nrows = 5;
else
    nrows = 4;
end
mdl_c = [27,158,119;...
         217,95,2;...
         117,112,179;...
         231,41,138;...
         102,166,30;...
         230,171,2;...
         166,118,29;...
         102,102,102] ./ 256;
mdl_clite = [102,194,165;...
             252,141,98;...
             141,160,203;...
             231,138,195;...
             166,216,84;...
             255,217,47;...
             229,196,148;...
             179,179,179] ./ 256;


% reshape predictions
y = mean(reshape(r.y,twin*ops.fs,[])');
c = mean(reshape(r.c,twin*ops.fs,[])');
cc = ops.mean_contrast ./ c;
cnt = 0;
for i = 1:size(r.pred,1)
    p(i,:) = mean(reshape(r.pred(i,:),twin*ops.fs,[])');
    if any(i == wI)
        cnt = cnt + 1;
        w(cnt,:) = mean(reshape(r.w(w_mdls == wI,:),twin*ops.fs,[])');
    end
end
mdls = ops.mdls;

% make contrast indices
c1 = nan(size(c)); c1(c==min(c)) = 1;
c2 = nan(size(c)); c2(c==max(c)) = 1;
time = (0:length(c)-1) ./ ops.fs;


% plot waveform
ax1 = subplot(nrows,ncols,3); hold on;
plot((1:length(u.waveform))/30e3*1e3,u.waveform*0.195,'k','linewidth',1);
plot([.2 .2],[-20 -40],'k','linewidth',2);
text(.25,-30,'20 \muV');
plot([.2 .7],[-45 -45],'k','linewidth',2);
text(.45,-50,'500 \mus','HorizontalAlignment','center');
colorbar;
title(sprintf('%s',res.cellID),'interpreter','none');
axis tight; axis square; plotPrefs;

%  % strf
%  map = [140,81,10;...
%         216,179,101;...
%         246,232,195;...
%         256,256,256;...
%         199,234,229;...
%         90,180,172;...
%         1,102,94] ./ 256;
%  ax2 = subplot(nrows,ncols,6);
%  strf = res.strf;
%  clims = linspace(-max(abs(strf(:))),max(abs(strf(:))),1000);
%  % cmap = zeroCMap(clims,0,map(1,:),map(4,:),map(end,:));
%  % cmap = zeroCMap(clims,0,[0.677, 0.492, 0.093],[1 1 1],[0.217,
%  % 0.525, 0.910]);
%  cmap = zeroCMap(clims,0);
%  imagesc(ops.t,ops.f/1000,strf);
%  h = colorbar; colormap(ax2,cmap);
%  caxis([clims(1) clims(end)]);
%  xlabel('Time (s)'); ylabel('Freq. (kHz)');
%  title('STRF');
%  plotPrefs; axis square;

% strf
ax2 = subplot(nrows,ncols,6);
strf = res.strf;
clims = linspace(-max(abs(strf(:))),max(abs(strf(:))),1000);
cmap = flipud(cbrewer2('PiYG'));
imagesc(ops.t,ops.f/1000,strf);
h = colorbar; colormap(ax2,cmap);
caxis([clims(1) clims(end)]);
xlabel('Time (s)'); ylabel('Freq. (kHz)');
title('STRF');
plotPrefs; axis square;

% nl per contrast
ax3 = subplot(nrows,ncols,9); hold on;
cols = {'b','r'};
cols_lite = {[.5 .5 1], [1 .5 .5]};
for i = 1:2
    xx = res.gc_ln.x{i}; yy = res.gc_ln.y{i};
    xf = linspace(min(xx),max(xx),100);
    scatter(xx,yy,5,cols{i},'.');
    plot(xf,ops.ln_model(res.gc_ln.ahat(i,:),xf),cols{i},'linewidth',1);
end
colorbar;
ylabel('spks/s');
xlabel('$X \cdot \beta_{strf}$','interpreter','latex');
axis tight; title('Nonlinearity');
axis square; plotPrefs;

ax1.Position(3:4) = ax2.Position(3:4);
ax3.Position(3:4) = ax2.Position(3:4);


% spike raster
ax4 = subplot(nrows,ncols,1:2); hold on;
edges = -.1:.005:twin+.1;
plot(time,c1*length(s.block.stimOn)/(twin/ops.blockLength)+2,...
     'color',cols{1},'linewidth',2);
plot(time,c2*length(s.block.stimOn)/(twin/ops.blockLength)+2,...
     'color',cols{2},'linewidth',2);
[~,raster,trials] = makePSTH(u.spikes,s.block.stimOn(1:(twin/ops.blockLength):end),edges);
[~,sortI] = sort(ops.order_r(1:(twin/ops.blockLength):end,2));
[~,sortTrials] = ismember(trials,sortI);
scatter(raster,sortTrials,5,'k.')
axis tight; xlim([0 twin]);
ylabel('Trial'); plotPrefs;
title('Raster'); 

% mean fr and predictions
smoothing = 2; clear ph;
ax5 = subplot(nrows,ncols,4:5); hold on;
plot(time,c1*max(p(:)),...
     'color',cols{1},'linewidth',2);
plot(time,c2*max(p(:)),...
     'color',cols{2},'linewidth',2);
ph(1) = plot(time,SmoothGaus(y,smoothing),'k-','LineWidth',1);
for i = 1:length(pI)
    ph(i+1) =  plot(time,SmoothGaus(p(pI(i),:),smoothing),'-','color',mdl_c(pI(i),:));
    text(10,max(p(:)) - range(p(:)).*(pI(i))*.075,...
         sprintf('r=%03.2f',r.corr_all(pI(i))),'Color',mdl_c(pI(i),:));
end
ylabel('Spike Count'); plotPrefs;
legend(ph,['spike rate',mdls(pI)],'location','northeastoutside','interpreter','none')
title('Prediction'); axis tight;

% gain estimate
ax6 = subplot(nrows,ncols,7:8); hold on; clear ph;
plot(time,c1*2,'color',cols{1},'linewidth',2);
plot(time,c2*2,'color',cols{2},'linewidth',2);
plot(time,cc,'k--');
plot(time,ones(size(cc)),'k:');
for i = 1:length(wI)
    if ~isfield(res,'basis')
        plot(time,w(i,:),':','color',mdl_c(wI(i),:));
        ph(i) = plot(time,SmoothGaus(w(i,:),smoothing),'-', ...
                     'linewidth',1,'color',mdl_c(wI(i),:));
    else
        ph(i) = plot(time,w(i,:),'-', ...
                     'linewidth',1,'color',mdl_c(wI(i),:));
    end
end
if length(wI) == 1
    %plot interaction predictor
    plot(time,cc.*normalize(SmoothGaus(p(1,:),smoothing))/10+1,'k-');
end
ylabel('w'); axis tight; ylim([0 2]);
plotPrefs; title('Gain Estimate and Contrast');
legend(ph,mdls(wI),'location','northeastoutside')

% w for each contrast
ax7 = subplot(nrows,ncols,10:11); hold on; clear ph;
lstyle = {':','--','-'};
ind = -10:75;
times = (ind) / ops.fs;
plot([times(1) times(end)],[1.5 1.5],'k--');
plot([times(1) times(end)],[.5 .5],'k--');
plot([times(1) times(end)],[1 1],'k:');
plot([0 0],[0 2],'k','linewidth',1)
for i = 1:length(wI)
    if ~isfield(res,'basis')
        plot(times,SmoothGaus(w(i,120+ind(1):120+ind(end)),2),...
             lstyle{wI(i)-w_mdls(1)+1},'color',cols{2},'linewidth',1);
        ph(i) = plot(times,SmoothGaus(w(mdl(i),240+ind(1):240+ind(end)),2), ...
                     lstyle{wI(i)-w_mdls(1)+1},'color',cols{1},'linewidth',1);
    else
        plot(times,w(i,120+ind(1):120+ind(end)),...
             lstyle{wI(i)-w_mdls(1)+1},'color',cols{2},'linewidth',1);
        ph(i) = plot(times,w(i,240+ind(1):240+ind(end)), ...
                     lstyle{wI(i)-w_mdls(1)+1},'color',cols{1},'linewidth',1);
    end
    
    % plot fit if only one model
    if length(wI) == 1
        x = times(find(times==0):end); xf = linspace(0,max(times),100);
        y = w(1,120:120+ind(end));
        [prms,model,tau] = fitExpGrid(x,y);
        plot(xf,model(prms,xf),'color',cols{2},'linewidth',2);
        y = w(1,240:240+ind(end));
        [prms,model,tau] = fitExpGrid(x,y);
        plot(xf,model(prms,xf),'color',cols{1},'linewidth',2);
    end
    
end
legend(ph,mdls(wI),'location','northeastoutside');
axis tight; ylim([0 2]);
xlabel('Time (s)'); ylabel('w');

if isfield(r,'b2_kernel')
    tim = (0:size(r.b2_kernel{1},1)-1)/ops.fs;
    for i = 1:length(wI)
        subplot(nrows,ncols,13+(i+1)); hold on;
        if size(r.b3_kernel{w_mdls==wI(i)},2) == 1
            plot(tim,r.b3_kernel{w_mdls==wI(i)},'k--');
        else
            plot(tim,r.b3_kernel{w_mdls==wI(i)}(:,1),'b--');
            plot(tim,r.b3_kernel{w_mdls==wI(i)}(:,2),'r--');
        end
        if size(r.b2_kernel{w_mdls==wI(i)},2) == 1
            plot(tim,r.b2_kernel{w_mdls==wI(i)},'k');
        else
            plot(tim,r.b2_kernel{w_mdls==wI(i)}(:,1),'b');
            plot(tim,r.b2_kernel{w_mdls==wI(i)}(:,2),'r');
        end
        if i == 1
            legend('Contrast Kernel','Gain Kernel');
        end
        xlabel('Time (s)'); ylabel('Kernel Weight');
        title(mdls{wI(i)});
    end
end

ax4.Position(3) = ax7.Position(3);
ax5.Position(3) = ax7.Position(3);
ax6.Position(3) = ax7.Position(3);
ax7.Position(3) = ax7.Position(3)*.5;

if length(wI) == 1
    % plot the stimulus
    ax8 = subplot(nrows,ncols,13:14);
    X = s.block.stimInfo.DB(:,1:length(time));
    %clims = linspace(min(X(:)),max(X(:)),1000);
    %cmap = zeroCMap(clims,ops.stimInfo.mu);
    imagesc(time,ops.f/1000,X); colorbar;
    %h = colorbar; colormap(ax8,cmap);
    %caxis([clims(1) clims(end)]);
    xlabel('Time (s)'); ylabel('Freq. (kHz)');
    plotPrefs;
    ax8.Position(3) = ax4.Position(3);
end
    

% weights
ax = subplot(nrows,ncols,12); hold on; clear xx;
for i = length(wI):-1:1
    x1 = [1 10 20 30+2*size(r(1).beta2{1},1)];
    w_i = w_mdls==wI(i);
    xx{1} = x1(1);
    xx{2} = x1(2);
    xx{3} = x1(3):x1(3)+numel(r.beta2{w_i})-1;
    xx{4} = x1(4):x1(4)+numel(r.beta3{w_i})-1;
    plot(xx{1},r.beta0{w_i},'-','color',mdl_c(wI(i),:),'marker','.');
    plot(xx{2},r.beta1{w_i},'-','color',mdl_c(wI(i),:),'marker','.');
    plot(xx{3},r.beta2{w_i},'-','color',mdl_c(wI(i),:),'marker','.');
    plot(xx{4},r.beta3{w_i},'-','color',mdl_c(wI(i),:),'marker','.');
end
set(gca,'xtick',[1 10 20 30+80])
set(gca,'xticklabels',{'b0','b1','b2','b3'});
%legend('b0','b1','b2','b3','location','northeastoutside');
title('Predictor Weights'); plotPrefs;
