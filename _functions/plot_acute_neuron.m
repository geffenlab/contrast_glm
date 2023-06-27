function plot_acute_neuron(filename)

load(filename)

ncols = 3;
nrows = 3;
blockLen = ops.stimInfo.baseNoiseD;
twin = blockLen * 4; % multiple of blockLen

% spectrogram supplied?
if exist('spec','var')
    spec_r = repmat(spec',nreps,1) - mean(spec(:));
    order_r = repmat(order,nreps,1);
    x1 = reshape(mean(spec_r,2),blockLen*ops.fs,[])';
    clims = linspace(-max(abs(x1(:))),max(abs(x1(:))),1000);
    cmap = zeroCMap(clims,0);
    % stimulus and repeats
    ax1 = subplot(nrows,ncols,1); hold on;
    [~,sortI] = sortrows(ops.order_r);
    imagesc(x1(sortI,:));
    h = colorbar; colormap(ax1,cmap);
    caxis([clims(1) clims(end)]);
    title('Stimulus Patterns','interpreter','none');
    axis tight; plotPrefs;
    ylabel('Trial');
end

% spike raster
subplot(nrows,ncols,4);
edges = 0:.005:3;
[~,raster,trials] = makePSTH(u.spikes,s.block.stimOn,edges);
[~,sortTrials] = ismember(trials,sortI);
scatter(raster,sortTrials,10,'k.')
xlabel('Time (samples)');
ylabel('Trial'); plotPrefs;
title('Raster');


% reshape by trial
cc = ops.mean_contrast ./ res.c;
p0 = reshape(res.strf_fit.pred,twin*ops.fs,[])';
p1 = reshape(exp(res.strf_fit.pred+res.strf_fit.coeffs(1)),twin*ops.fs,[])';
p2 = reshape(res.gain_fit.pred,twin*ops.fs,[])';
p3 = reshape(res.gc_ln.pred,twin*ops.fs,[])';
y1 = reshape(res.y(ops.include),twin*ops.fs,[])';
c1 = mean(reshape(res.c(ops.include),twin*ops.fs,[])');
cc1 = mean(reshape(cc(ops.include),twin*ops.fs,[])');
w1 = reshape(res.w,twin*ops.fs,[])';


% strf
ax2 = subplot(nrows,ncols,2);
strf = res.strf;
clims = linspace(-max(abs(strf(:))),max(abs(strf(:))),1000);
cmap = zeroCMap(clims,0);
imagesc(ops.t,ops.f/1000,strf);
h = colorbar; colormap(ax2,cmap);
caxis([clims(1) clims(end)]);
xlabel('Time (s)'); ylabel('Freq. (kHz)');
title(sprintf('%s\nSTRF',res.cellID),'interpreter','none');
plotPrefs;

% nl per contrast
subplot(nrows,ncols,3); hold on;
cols = {'b','r'};
cols_lite = {[.5 .5 1], [1 .5 .5]};
for i = 1:2
    xx = res.gc_ln.x{i}; yy = res.gc_ln.y{i};
    xf = linspace(min(xx),max(xx),100);
    scatter(xx,yy,cols{i},'.');
    plot(xf,ops.ln_model(res.gc_ln.ahat(i,:),xf),cols{i},'linewidth',1);
end
ylabel('spks/s');
xlabel('$X \cdot \beta_{strf}$','interpreter','latex');
axis tight; title('Nonlinearity');
plotPrefs;


% mean fr and predictions
subplot(nrows,ncols,5:6); hold on;
yyaxis right;
plot(c1,'r','LineWidth',1);
ylabel('Contrast');
yyaxis left;
ph(1) = plot(SmoothGaus(mean(y1),1),'k-','LineWidth',1);
ph(2) = plot(SmoothGaus(mean(p1),1),'-','color',[.5 .5 .5]);
ph(3) = plot(SmoothGaus(mean(p3),1),'c-');
ph(4) = plot(SmoothGaus(mean(p2),1),'m-');
ylabel('Spike Count'); plotPrefs;
ax = gca; ax.YAxis(1).Color = 'k'; 
ax.YAxis(2).Color = 'r';
legend(ph,'spike rate','strf model','gain LN','gain GLM','location','ne')
title('Prediction'); axis tight;

% gain estimate
subplot(nrows,ncols,8:9); hold on;
yyaxis left
plot(cc1,'--','Color',[.5 .5 .5]);
plot(mean(w1),'k:');
plot(SmoothGaus(mean(w1),5),'k','linewidth',1);
ylabel('w'); ylim([0 2]);
yyaxis right
plot(c1,'r','linewidth',1);
ylabel('contrast');
axis tight; plotPrefs;
ax = gca; ax.YAxis(1).Color = 'k'; 
ax.YAxis(2).Color = 'r';
title('Gain Estimate and Contrast');
xlabel('Time (samples)'); 
colors = {'k',[.4 .4 .4],'m'};

% weights
ax = subplot(nrows,ncols,[7]); hold on; clear xx;
x1 = [1 10 20 30+numel(res.beta2)];
xx{1} = x1(1);
xx{2} = x1(2);
xx{3} = x1(3):x1(3)+numel(res.beta2)-1;
xx{4} = x1(4):x1(4)+numel(res.beta3)-1;
if isfield(res.gain_fit,'stats')
    errorbar(xx{1},res.beta0,res.gain_fit.stats.se(1),...
             'linewidth',.5,'color',[.8 .8 .8]);
    errorbar(xx{2},res.beta1,res.gain_fit.stats.se(2),...
             'linewidth',.5,'color',[.8 .8 .8]);
    errorbar(xx{3},...
             res.beta2,...
             res.gain_fit.stats.se(3:3+numel(res.beta2)-1),...
             'linewidth',.5,'color',[.8 .8 .8]);
    errorbar(xx{4},...
             res.beta3,...
             res.gain_fit.stats.se(3+numel(res.beta2):end),...
             'linewidth',.5,'color',[.8 .8 .8]);
end
if numel(res.beta2) == 2*ops.gain_lags
    plot(xx{1},res.beta0,'-','color',colors{1},'marker','.');
    plot(xx{2},res.beta1,'-','color',colors{2},'marker','.');
    plot(xx{3}(1:ops.gain_lags),res.beta2(1:ops.gain_lags),...
         '-','color',cols{1},'marker','.');
    plot(xx{3}(ops.gain_lags+1:end),res.beta2(ops.gain_lags+1:end),...
         '-','color',cols{2},'marker','.');
else
    plot(xx{1},res.beta0,'-','marker','.');
    plot(xx{2},res.beta1,'-','marker','.');
    plot(xx{3},res.beta2,'-','marker','.');
end
if numel(res.beta3) == 2*ops.gain_lags
    plot(xx{4}(1:ops.gain_lags),res.beta3(1:ops.gain_lags),...
         '-','color',cols_lite{1},'marker','.');
    plot(xx{4}(ops.gain_lags+1:end),res.beta3(ops.gain_lags+1:end),...
         '-','color',cols_lite{2},'marker','.');
else
    plot(xx{4},res.beta3,'-','marker','.','color',colors{3});
end
%symlog(ax,'y',-1.5); set(ax,['y','MinorGrid'],'off');
set(gca,'xtick',[1 10 20 30+numel(res.beta2)])
set(gca,'xticklabels',{'b0','b1','b2','b3'});
%legend('b0','b1','b2','b3','location','northeastoutside');
title('Predictor Weights'); plotPrefs;