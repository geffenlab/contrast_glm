function f = plot_sim_neuron(fn);

load(fn)

% get contrast
c = neuron.ops.contrast;

% trial samples
nts = (ops.blockLength*ops.fs*2 + ops.pad*ops.fs);

% number of trials to plot
nt = 2;

f = figure; set(f,'Position',[0 0 1200 400]);
subplot(2,5,1); hold on;
imagesc(ops.t,ops.f/1000,neuron.ops.beta); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('STRF (beta)'); axis tight; plotPrefs;

subplot(2,5,2); hold on;
scatter(neuron.h(c==max(ops.sd)),neuron.y(c==max(ops.sd)),'r.');
scatter(neuron.h(c==min(ops.sd)),neuron.y(c==min(ops.sd)),'b.');
xlabel('Linear Drive');
ylabel('Firing Rate');
title('Nonlinearity'); plotPrefs;

subplot(2,5,3); hold on;
imagesc(reshape(neuron.h,nts,[])'); colorbar;
plot(repmat(ops.blockLength*ops.fs,2,1),[0 ops.ntrials*ops.exemplars],'r','linewidth',2);
xlabel('Time (step)');
ylabel('Trial');
title(sprintf('Linear Drive\nc_{low}=%02.1f,c_{hi}=%02.1f',ops.sd(1),ops.sd(2)));
axis tight;
plotPrefs;

subplot(2,5,4); hold on;
imagesc(reshape(neuron.y,nts,[])'); colorbar;
plot(repmat(ops.blockLength*ops.fs,2,1),[0 ops.ntrials*ops.exemplars],'r','linewidth',2);
plot(mean(reshape(neuron.y,nts,[])'),'w','linewidth',1)
xlabel('Time (step)');
ylabel('Trial');
title('Poisson Counts');
axis tight; plotPrefs;

subplot(2,5,5); hold on;
yyaxis right;
plot(mean(reshape(neuron.ops.contrast,nts*nt,[])'),'color','r');
ylabel('Contrast'); 
ylim(ylim + [-1 1]);
yyaxis left;
plot(SmoothGaus(mean(reshape(neuron.y,nts*nt,[])'),2),'k','LineWidth',1);
ylabel('Spike Count');
xlabel('Time (step)');
title('Spike Rate');
xlim([0 nts*nt]);
plotPrefs;

subplot(2,5,10); hold on;
inst_gain = mean(reshape(neuron.g,nts*nt,[])',1);
plot(inst_gain,'k','LineWidth',1);
xlabel('Time (step)');
ylabel('Gain');
title('Gain');
plotPrefs;

if exist('res','var') && ~isempty(res)
    
    b1 = res.beta1;
    b2 = res.beta2;
    
    subplot(2,5,5); hold on;
    plot(SmoothGaus(mean(reshape(res.gain_fit.pred,nts*nt,[])'),2),'m','LineWidth',1)
    
    subplot(2,5,6); hold on;
    imagesc(ops.t,ops.f/1000,res.strf); colorbar;
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    title('Fitted STRF'); axis tight; plotPrefs;
    
    subplot(2,5,7:8); hold on;
    if isfield(res.gain_fit,'stats')
        errorbar(1,res.gain_fit.stats.beta(1),res.gain_fit.stats.se(1),'marker','.');
        errorbar(3,res.gain_fit.stats.beta(2),res.gain_fit.stats.se(2),'marker','.');
        errorbar(5:5+length(res.beta2)-1,...
                 res.gain_fit.stats.beta(3:3+length(res.beta2)-1),...
                 res.gain_fit.stats.se(3:3+length(res.beta2)-1),'marker','.');
        errorbar(7+length(res.beta2):7+length(res.beta2)+length(res.beta3)-1,...
                 res.gain_fit.stats.beta(3+length(res.beta2):end),...
                 res.gain_fit.stats.se(3+length(res.beta2):end),'marker','.');
    else
        plot(1,res.beta0,'marker','.');
        plot(3,res.beta1,'marker','.');
        plot(5:5+length(res.beta2)-1,res.beta2,'marker','.');
        plot(7+length(res.beta2):7+length(res.beta2)+length(res.beta3)-1,...
             res.beta3,'marker','.');
    end
    legend('b0','b1','b2','b3','location','se');
    title('Predictor Weights'); plotPrefs;
    
    if isfield(res,'w1')
        % "combined" w
        subplot(2,5,9); hold on;
        plot(1:160,neuron.g(1:160),'k')
        plot(1:160,[res.w1(1:80) res.w2(81:160)],'--','color',[.5 .5 .5])
        plot(1:80,res.w1(1:80),'.b');
        plot(81:160,res.w2(81:160),'.r');
        xlabel('Time (step)');
        ylabel('Gain');
        legend('Gain','Dummy-Coded Estimates');
    elseif isfield(res,'basis')
        subplot(2,5,9); hold on;
        nb = sum(ops.spline_basis)-2;
        plot((1:ops.gain_lags) / ops.fs,res.basis*res.beta2(1:nb),'b','linewidth',1)
        plot((1:ops.gain_lags) / ops.fs,res.basis*res.beta2(nb+1:end),'r','linewidth',1);
        ylabel('Filter Weight'); xlabel('Time (s)');
        title('Contrast Kernels');
        legend('Low Contrast','High Contrast');
        
    end

    
    subplot(2,5,10); hold on;
    if isfield(res,'wm')
        w = reshape(res.wm',1,numel(res.wm));
    elseif isfield(res,'w')
        w = res.w;
    end
    wmt = reshape(w,nts*nt,[])';
    plot(mean(wmt,1),'r','LineWidth',1);        
    if ~isfield(res,'ind_model')
        legend('True','Estimate');
    else
        plot(res.ind_model.w_norm(1:nts*nt),'m','LineWidth',1);
        legend('True','Full Model','Ind. Model');
    end
    axis tight; plotPrefs;
        
else
    
    legend('True');

end






