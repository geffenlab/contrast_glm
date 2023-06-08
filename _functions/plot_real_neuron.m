function f = plot_real_neuron(filename)

%% setup
load(filename)

% restate all results in original timescale (ommitted data isnan)
pp1 = nan(size(ops.include));
pp1(ops.include) = res.strf_fit.pred;
pp2 = pp1;
pp2(ops.include) = res.gain_fit.pred;

% get w estimate
ww = res.w;

% get spikes
yy = res.y;

% get contrast
c = res.contrast;

% reshape data according to trials (sample up to min length)
trial_starts = find(diff([0; res.trialI(:,1)]));
mntl = min(res.stimLength);

for i = 1:length(trial_starts)
    I = trial_starts(i):(trial_starts(i)+mntl);
    w(i,:) = ww(I);
    y(i,:) = yy(I);
    cc(i,:) = c(I);
    p1(i,:) = pp1(I);
    p2(i,:) = pp2(I);
end
[~,sortI] = sort(d.trialType(:,3));



%% plotting
f = figure(9124); clf;

% plot strf
subplot(2,6,1); hold on;
imagesc(ops.t,ops.f/1000,reshape(res.strf,length(ops.f),[])); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('STRF'); axis tight;

% plot nonlinearities
subplot(2,6,2); hold on;
scatter(pp1(c==max(ops.contrast)),yy(c==max(ops.contrast)),'r.');
scatter(pp1(c==min(ops.contrast)),yy(c==min(ops.contrast)),'b.');
xlabel('Linear Drive');
ylabel('Firing Rate');
title('Nonlinearity')

% plot linear drive
subplot(2,6,3:4); hold on;
imagesc(p1(sortI,:)); colorbar;
plot(repmat(3*ops.fs,2,1),[1 size(w,1)],'r','linewidth',2)
xlabel('Time (step)');
ylabel('Trial');
title(sprintf('Linear Drive\nc_{low}=%02.1f,c_{hi}=%02.1f',ops.contrast(1),ops.contrast(2)));
axis tight;

% plot spike raster
subplot(2,6,5:6); hold on;
imagesc(y(sortI,:)); colorbar;
plot(repmat(3*ops.fs,2,1),[1 size(w,1)],'r','linewidth',2)
plot(mean(y),'w','linewidth',1)
xlabel('Time (step)');
ylabel('Trial');
title('Spike Count');
axis tight;

% plot coefficients
subplot(2,6,7:8); hold on
errorbar(1,res.gain_fit.stats.beta(1),res.gain_fit.stats.se(1),'marker','.');
errorbar(3,res.gain_fit.stats.beta(2),res.gain_fit.stats.se(2),'marker','.');
errorbar(5:5+length(res.beta3)-1,...
    res.gain_fit.stats.beta(3:3+length(res.beta3)-1),...
    res.gain_fit.stats.se(3:3+length(res.beta3)-1),'marker','.');
errorbar(7+length(res.beta3):7+length(res.beta3)*2-1,...
    res.gain_fit.stats.beta(3+length(res.beta3):end),...
    res.gain_fit.stats.se(3+length(res.beta3):end),'marker','.');
legend('b0','b1','b2','b3','location','se'); title('Coefficients');


% plot fr, predictions
subplot(2,6,9:10); hold on;
yyaxis left;
plot(mean(y,'omitnan') * ops.fs,'k-','linewidth',1);
plot(mean(exp(p1 + res.strf_fit.coeffs(1)),'omitnan') * ops.fs,...
     '-','color',[.5 .5 .5]);
plot(mean(p2,'omitnan') * ops.fs,'m-');
ylabel('spks/s');
yyaxis right;
plot(mean(cc,'omitnan'),'r','LineWidth',1);
ylabel('Contrast');
xlabel('Time (step)');
legend('data','strf','gain glm','contrast');
ax = gca; ax.YAxis(1).Color = 'k'; 
ax.YAxis(2).Color = 'r'; 
title('Predictions'); axis tight;

% plot w, contrast
subplot(2,6,11:12); hold on;
yyaxis left
plot([1 size(w,2)],[1 1],'--','color',[.5 .5 .5]);
plot(mean(ops.mean_contrast./cc,'omitnan'),'k--');
plot(mean(w,'omitnan'),'k','LineWidth',1.5)
ylabel('w');
yyaxis right
plot(mean(log(cc)),'r-','LineWidth',1);
ylabel('Contrast');
ax = gca; ax.YAxis(1).Color = 'k'; 
ax.YAxis(2).Color = 'r'; 
xlabel('Time (step)');
title('Gain'); axis tight;










