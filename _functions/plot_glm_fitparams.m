function plot_glm_fitparams(d,include,ops,errtype)

if ~exist('errtype','var') | isempty(errtype);
    errtype = 'prctile';
end

nrows = 4;
ncols = 4;

% number of cells per condition
subplot(nrows,ncols,1)
p = pie([sum(d.contrastI==0&include),sum(d.contrastI==1&include)],...
        {'High-Low','Low-High'});
pcs = [sum(d.contrastI==0),sum(d.contrastI==1)] ./ ...
      length(d.contrastI) * 100;
p(1).FaceColor = 'b';
p(2).String = sprintf('%s\n%3.2f%%\n(n = %d)',...
                      p(2).String,pcs(1),sum(d.contrastI==0&include));
p(4).String = sprintf('%s\n%3.2f%%\n(n = %d)',...
                      p(4).String,pcs(2),sum(d.contrastI==1&include));
p(3).FaceColor = 'r';

% prediction strength (correlation)
subplot(nrows,ncols,2); hold on;
edges = linspace(0,1,50);
centers = edges(1:end-1) + mean(diff(edges));
histogram(d.p1e_corr(d.contrastI == 0 & include),edges,...
          'FaceColor','b','FaceAlpha',.1);
histogram(d.p1e_corr(d.contrastI == 1 & include),edges,...
          'FaceColor','r','FaceAlpha',.1);
histogram(d.p2_corr(d.contrastI == 0 & include),edges,...
          'FaceColor','b');
histogram(d.p2_corr(d.contrastI == 1 & include),edges,...
          'FaceColor','r');
xlabel('r'); ylabel('Cells'); xlim([0 1]); plotPrefs;
legend('H-L (strf)','L-H (strf)','H-L (gain)','H-L (gain)',...
       'location','nw');
title('Predictive Power'); plotPrefs;

% plot mean spike rate
subplot(nrows,ncols,5:6); hold on; clear p;
t1 = 1:3*ops.fs;
t2 = 3*ops.fs:size(d.y_new,2);
plot([120 120],[0 1],'k--','LineWidth',1.5);
p(1) = plot(t1,mean(d.p2_new(include&d.contrastI==0,t1)),'r--','LineWidth',1);
plot(t2,mean(d.p2_new(include&d.contrastI==0,t2)),'b--','LineWidth',1);
plot(t1,mean(d.p2_new(include&d.contrastI==1,t1)),'b--','LineWidth',1);
plot(t2,mean(d.p2_new(include&d.contrastI==1,t2)),'r--','LineWidth',1);
p(2) = plot(t1,mean(d.p1e_new(include&d.contrastI==0,t1)),'r:','LineWidth',1);
plot(t2,mean(d.p1e_new(include&d.contrastI==0,t2)),'b:','LineWidth',1);
plot(t1,mean(d.p1e_new(include&d.contrastI==1,t1)),'b:','LineWidth',1);
plot(t2,mean(d.p1e_new(include&d.contrastI==1,t2)),'r:','LineWidth',1);
[~,p(3)] = patchErrorBars(t1,d.y_new(include&d.contrastI==0,t1),'r',errtype);
patchErrorBars(t2,d.y_new(include&d.contrastI==0,t2),'b',errtype);
patchErrorBars(t1,d.y_new(include&d.contrastI==1,t1),'b',errtype);
patchErrorBars(t2,d.y_new(include&d.contrastI==1,t2),'r',errtype);
axis tight; ylabel('Spike Count'); xlabel('Time (step)');
legend(p,{'Gain Model Fit','STRF Model Fit','Observed Spikes'},...
       'location','nw');
title('Spike Rate and Predictions'); plotPrefs;

% plot mean spike rate
subplot(nrows,ncols,9); hold on; clear p;
t1 = 1:3*ops.fs;
t2 = 3*ops.fs:size(d.y_new,2);
p(1) = plot(t1,mean(d.p2_new(include&d.contrastI==0,t1)),'r--','LineWidth',1);
plot(t2,mean(d.p2_new(include&d.contrastI==0,t2)),'b--','LineWidth',1);
plot(t1,mean(d.p2_new(include&d.contrastI==1,t1)),'b--','LineWidth',1);
plot(t2,mean(d.p2_new(include&d.contrastI==1,t2)),'r--','LineWidth',1);
p(2) = plot(t1,mean(d.p1e_new(include&d.contrastI==0,t1)),'r:','LineWidth',1);
plot(t2,mean(d.p1e_new(include&d.contrastI==0,t2)),'b:','LineWidth',1);
plot(t1,mean(d.p1e_new(include&d.contrastI==1,t1)),'b:','LineWidth',1);
plot(t2,mean(d.p1e_new(include&d.contrastI==1,t2)),'r:','LineWidth',1);
[~,p(3)] = patchErrorBars(t1,d.y_new(include&d.contrastI==0,t1),'r',errtype);
patchErrorBars(t2,d.y_new(include&d.contrastI==0,t2),'b',errtype);
patchErrorBars(t1,d.y_new(include&d.contrastI==1,t1),'b',errtype);
patchErrorBars(t2,d.y_new(include&d.contrastI==1,t2),'r',errtype);
axis tight; xlim([0 40]); ylim([0 1]); plotPrefs;

% plot mean spike rate
subplot(nrows,ncols,10); hold on; clear p;
t1 = 1:3*ops.fs;
t2 = 3*ops.fs:size(d.y_new,2);
plot([120 120],[0 1],'k--','LineWidth',1.5);
p(1) = plot(t1,mean(d.p2_new(include&d.contrastI==0,t1)),'r--','LineWidth',1);
plot(t2,mean(d.p2_new(include&d.contrastI==0,t2)),'b--','LineWidth',1);
plot(t1,mean(d.p2_new(include&d.contrastI==1,t1)),'b--','LineWidth',1);
plot(t2,mean(d.p2_new(include&d.contrastI==1,t2)),'r--','LineWidth',1);
p(2) = plot(t1,mean(d.p1e_new(include&d.contrastI==0,t1)),'r:','LineWidth',1);
plot(t2,mean(d.p1e_new(include&d.contrastI==0,t2)),'b:','LineWidth',1);
plot(t1,mean(d.p1e_new(include&d.contrastI==1,t1)),'b:','LineWidth',1);
plot(t2,mean(d.p1e_new(include&d.contrastI==1,t2)),'r:','LineWidth',1);
[~,p(3)] = patchErrorBars(t1,d.y_new(include&d.contrastI==0,t1),'r',errtype);
patchErrorBars(t2,d.y_new(include&d.contrastI==0,t2),'b',errtype);
patchErrorBars(t1,d.y_new(include&d.contrastI==1,t1),'b',errtype);
patchErrorBars(t2,d.y_new(include&d.contrastI==1,t2),'r',errtype);
axis tight; xlim([120 160]); ylim([0 1]); plotPrefs;

% plot w traces
subplot(nrows,ncols,13:14); hold on;
patchErrorBars(t1,d.w_new(include&d.contrastI==0,t1),'r',errtype);
patchErrorBars(t2,d.w_new(include&d.contrastI==0,t2),'b',errtype);
patchErrorBars(t1,d.w_new(include&d.contrastI==1,t1),'b',errtype);
patchErrorBars(t2,d.w_new(include&d.contrastI==1,t2),'r',errtype);
axis tight; plotPrefs; title('Gain Estimate'); ylabel('w');
plot([120 120],ylim,'k--','LineWidth',1.5);




% plot average centered strfs
subplot(nrows,ncols,3);
octs = [-ops.bins:ops.bins]/8;
imagesc(ops.t,octs,nanmean(d.cstrf(:,:,include&d.contrastI==0),3));
colorbar; xlabel('Lag (s)'); ylabel('Octaves from CF');
title('High-Low Contrast STRF'); plotPrefs;
subplot(nrows,ncols,4);
imagesc(ops.t,octs,nanmean(d.cstrf(:,:,include&d.contrastI==1),3));
colorbar; xlabel('Lag (s)'); ylabel('Octaves from CF');
title('Low-High Contrast STRF'); plotPrefs;

% plot b0
subplot(nrows,ncols,7); hold on;
edges = linspace(min(d.b0(include)),max(d.b0(include)),50);
histogram(d.b0(include&d.contrastI==0),edges,'FaceColor','b');
histogram(d.b0(include&d.contrastI==1),edges,'FaceColor','r');
xlabel('${\beta_0}$','Interpreter','latex');
title('Baseline'); plotPrefs;


% plot b1
subplot(nrows,ncols,8); hold on;
edges = linspace(min(d.b1(include)),max(d.b1(include)),50);
histogram(d.b1(include&d.contrastI==0),edges,'FaceColor','b');
histogram(d.b1(include&d.contrastI==1),edges,'FaceColor','r');
xlabel('${\beta_1}$','Interpreter','latex');
title('Stimulus Scaling'); plotPrefs;


% plot b2
subplot(nrows,ncols,11:12); hold on;
patchErrorBars(1:ops.gain_lags,d.b2(include&d.contrastI==0,:),'b',errtype);
patchErrorBars(1:ops.gain_lags,d.b2(include&d.contrastI==1,:),'r',errtype);
ylabel('${\beta_2}$','Interpreter','latex'); xlim([1 ops.gain_lags]);
title('Contrast x Stimulus Interaction'); plotPrefs;


% plot b3
subplot(nrows,ncols,15:16); hold on;
patchErrorBars(1:ops.gain_lags,d.b3(include&d.contrastI==0,:),'b',errtype);
patchErrorBars(1:ops.gain_lags,d.b3(include&d.contrastI==1,:),'r',errtype);
ylabel('${\beta_3}$','Interpreter','latex'); xlabel('Lag (step)');
title('Contrast Scaling'); plotPrefs; xlim([1 ops.gain_lags]);