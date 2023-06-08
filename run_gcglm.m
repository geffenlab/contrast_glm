%clear all; close all;
addpath(genpath('./_functions/'));

% to run, download the dataset from Dryad and specify the path to the GLM results here:
%resPath = '~/data/gain_opto/glm_res';
dc = [0 1 2];
models = {'ln-static','ln-gc','glm','glm-symm','glm-asymm(1)', ...
          'glm-asymm(2)'};
ops.mdls = models;
tag = 'spline-3-3_';
resfile = sprintf('./acute_res_%s2tau_lohse.mat',tag);

% setup figure
plot_on = false;
plot_visible = false;
sz = [900 900];
f1 = figure(1); set(f1,'visible',plot_visible,'Position',[0 0 sz]);


if ~exist(resfile,'file')
    % get file list for one model (all neurons)
    fileList = dir(fullfile(resPath,sprintf('acute_%sdc%d*.mat',tag,dc(1))));

    for i = 1:length(fileList)
        
        fprintf('%d/%d... ',i,length(fileList)); tic;
        
        for m = 1:length(dc)
            
            % get file list for this model
            fileList = dir(fullfile(resPath,sprintf('acute_%sdc%d*.mat',tag,dc(m))));

            load(fullfile(fileList(i).folder,fileList(i).name));
            
            % common variabiles
            if m == 1
                
                % centered strf with lag/bf
                ops.bins = 7;
                [cstrf,bf,lag] = centeredSTRF(res.strf,ops);
                
                % activity/strf
                r(i).cellID = res.cellID;
                r(i).sessionID = res.sessionID;
                r(i).y = res.y;
                r(i).c = res.c;
                r(i).scene = res.scene;
                r(i).strf = res.strf;
                r(i).cstrf = cstrf;
                r(i).bf = bf;
                r(i).lag = lag;
                r(i).nr = median(res.nr,1,'omitnan');
                
                % ln results
                r(i).static_ln_gain = res.static_ln.ahat(3);
                r(i).static_ln_mxslope = res.static_ln.maxslope;
                r(i).static_ln_meanslope = res.static_ln.maxslope;
                r(i).static_ln_beta = res.static_ln.beta;
                r(i).pred(1,:) = res.static_ln.pred;
                r(i).gc_ln_gain = res.gc_ln.ahat(:,3);
                r(i).gc_ln_mxslope = res.gc_ln.maxslope;
                r(i).gc_ln_meanslope = res.gc_ln.maxslope;
                r(i).gc_ln_beta = res.gc_ln.beta;
                r(i).pred(2,:) = res.gc_ln.pred;
                
                % static strf result (exponent)
                r(i).pred(3,:) = exp(res.strf_fit.pred + res.strf_fit.coeffs(1));
                
                % reshape STRF prediction to original order
                x = nan(size(ops.include));
                for j = 1:ops.cvfolds
                    x(ops.cv_sample ~= j) = res.strf_fit_cv(j).pred;
                end
                                
                % fit model from Lohse et al. (2020)                
                twin = 4 * ops.stimInfo.baseNoiseD;
                x1 = mean(reshape(x,twin*ops.fs,[])',1);
                y1 = mean(reshape(res.y(ops.include),twin*ops.fs,[])',1) * ops.fs;
                tau_max = 40; % max 40 gives close to 35% epoch (1 sec)
                tau_0 = 5; 
                ops.valid_conv = ops.blockLength * ops.fs;
                options = optimoptions('fmincon');
                options.Display = 'none';
                [tau_h,fval,exitflag,output] = fmincon(...
                    @(tau) erfTauNL(x1,y1,tau,res,ops),...
                    tau_0,[],[],[],[],1,tau_max,[],options);
                                
                % save results
                r(i).pred(4,:) = repmat(tauNL(x1,tau_h,res,ops),1, ...
                                        length(x)/length(x1)) / ops.fs;
                r(i).lohse_tau = tau_h;
                fprintf(' tau = %03.2f',tau_h);
                
            end
            
            % model results
            r(i).beta0{m} = res.beta0;
            r(i).beta1{m} = res.beta1;
            r(i).beta2{m} = res.beta2;
            r(i).beta3{m} = res.beta3;
            r(i).w(m,:) = res.w;
            r(i).pred(4+m,:) = res.gain_fit.pred;
            
            % reshape w (by 4 trials to get all transitions with padding)
            r(i).ws(m,:) = mean(reshape(res.w,ops.blockLength*4*ops.fs,[])',1);
            
            if isfield(res,'basis')
                nb = size(res.basis,2);
                if m == 1
                    r(i).b2_kernel{m} = res.basis*res.beta2;
                    r(i).b3_kernel{m} = res.basis*res.beta3;
                elseif m == 2
                    r(i).b2_kernel{m} = [res.basis*res.beta2(1:nb),...
                                        res.basis*res.beta2(nb+1:end)];
                    r(i).b3_kernel{m} = res.basis*res.beta3;
                elseif m == 3
                    r(i).b2_kernel{m} = [res.basis*res.beta2(1:nb),...
                                        res.basis*res.beta2(nb+1:end)];
                    r(i).b3_kernel{m} = [res.basis*res.beta3(1:nb),...
                                        res.basis*res.beta3(nb+1:end)];
                end
            end
            
        end
        
        % reshape y
        r(i).ys = reshape(r(i).y,ops.blockLength*ops.fs,[])';        
        
        % compute corr and r2 for each model for each contrast
        for ii = 1:size(r(i).pred,1)
            r(i).ps(ii,:,:) = reshape(r(i).pred(ii,:),ops.blockLength*ops.fs,[])';
            for jj = 1:2
                I = ops.order_r(:,1) == ops.contrast(jj);
                r(i).corr(ii,jj) = ...
                    corr(mean(squeeze(r(i).ys(I,:)))',...
                         mean(squeeze(r(i).ps(ii,I,:)))');
                r(i).r2(ii,jj) = r(i).corr(ii,jj).^2;
            end
        end
        
        % compute corr and r2 across contrasts
        yt = reshape(r(i).y,ops.blockLength*ops.fs*2,[])';
        for ii = 1:size(r(i).pred,1)
            pt = reshape(r(i).pred(ii,:),ops.blockLength*ops.fs*2,[])';
            r(i).corr_all(ii) = corr(mean(yt)',mean(pt)');
            r(i).r2_all(ii) = r(i).corr(ii).^2;
        end
        
        % compute excluding transition window
        I = repmat([false(1,ops.fs) true(1,ops.fs*(ops.blockLength-1))],1,2);
        for ii = 1:size(r(i).pred,1)
            pt = reshape(r(i).pred(ii,:),ops.blockLength*ops.fs*2,[])';
            r(i).corr_not(ii) = corr(mean(yt(:,I))',mean(pt(:,I))');
            r(i).r2_not(ii) = r(i).corr_not(ii).^2;
        end
        
        % cell stats
        r(i).mfr = u.mfr;
        r(i).trough_peak = u.trough_peak;
        r(i).peak_inflect = u.peak_inflect;
        r(i).FWHM = u.FWHM;
        
        % plot individual results
        if plot_on
            blks = 4; % number of blocks to repeat
            twin = blks * ops.stimInfo.baseNoiseD;
            ops.mdls = models;
            plot_acute_neuron_mdls(res,r(i),u,s,ops,twin);
            fn = sprintf('./_plots/%s.pdf',[tag r(i).cellID]);
            saveFigPDF(f1,sz,fn);
            clf(f1);
        end
        
        toc;
        
    end
    
    % model labels
    ops.mdls = models;
    ops.dc = dc;

    save(resfile,'r','ops');

else
    
    load(resfile);

end


%% post processing

% compute gain control metrics
ln_beta = cat(3,r.gc_ln_beta); ln_beta = squeeze(ln_beta(:,2,:));
%ln_gain = [r.gc_ln_gain];
gc(1,:) = ln_beta(2,:) - ln_beta(1,:);
glm_w = cat(3,r.ws);
glm_gain = glm_w(:,[116 234],:);
gc(2:4,:) = squeeze(glm_gain(:,2,:) - glm_gain(:,1,:));



%% inclusion criteria
nr = vertcat(r.nr)';
nr_crit = all(nr > 0 & nr < 50,1);
fr_crit = [r.mfr] > 1;
gc_crit = 0;
gc_metric = 4;
w_metric = gc_metric - 1;
cell_include = fr_crit & nr_crit;
gc_cut = prctile(gc(gc_metric,:),25); % use lowest 25% for adaptation??


%% fit gain timecourse
% gain time window for fitting
ind = -3:40;
cind = [240 120];
start = 0;
xtimes = (ind) / ops.fs; xf = linspace(xtimes(ind==start),xtimes(end),100);

% fit each neuron with an exponential for each contrast
true_tau = []; true_prms = [];
fprintf('Fitting exponentials... ');
for i = 1:size(gc,2)
    fprintf('%g ',i);
    for j = 1:2
        y = glm_w(w_metric,ind+cind(j),i);
        [true_prms(:,i,j),mdl,true_tau(i,j)] = fitExpGrid(xtimes(ind>=start),y(ind>=start));
    end
end
fprintf('\n\n');



        
%% plot stuff
c = {'b','r'};
c_rgb = [0 0 1; 1 0 0];
c_lite = [.6 .6 1; 1 .6 .6];
nrows = 3;
ncols = 3;
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
sz = [700 600];

%%%%%%%%%%%
%% plotting
f1 = figure(12); clf; set(f1,'Position',[0 0 sz]);

%% plot model performance
% whole time-course
ops.mdls = {'ln_static','ln-gc','glm-static','lohse','glm-symm','glm-asymm(1)','glm-asymm(2)'};
ax1 = subplot(nrows,ncols,1); hold on;
mdlI = [1 2 4 7];
plot([0 length(mdlI)+1],[0 0],'k');
corr_all = cat(1,r.corr_all);
plotDist(1:length(mdlI),corr_all(:,mdlI),num2cell(mdl_c(mdlI,:),2),...
         [],'prctile','median',[],'Marker','.',...
         'MarkerSize',5);
set(gca,'xtick',1:length(mdlI),'xticklabels',ops.mdls(mdlI));
xtickangle(45);
ylabel('r'); plotPrefs; axis tight; ylim([-1 1]);
title('Model Performance (all timepoints)');

% exclude adaptation period
ax2 = subplot(nrows,ncols,2); hold on;
plot([0 length(mdlI)+1],[0 0],'k');
corr_not = cat(1,r.corr_not);
plotDist(1:length(mdlI),corr_not(:,mdlI),num2cell(mdl_clite(mdlI,:),2),...
         [],'prctile','median',[],'Marker','.',...
         'MarkerSize',5);
set(gca,'xtick',1:length(mdlI),'xticklabels',ops.mdls(mdlI));
xtickangle(45);
ylabel('r'); plotPrefs; axis tight; ylim([-1 1]);
title('Model Performance (1s after switch)');


%x1 = [1 2 4 5 7 8 10 11 13 14];
%plotDist(x1,r2s,{'b','r','b','r','b','r','b','r','b','r'},...
%         [],'prctile','median',[],'Marker','.','MarkerSize',5);
%set(gca,'xtick',[1.5 4.5 7.5 10.5 13.5],'xticklabels',ops.mdls);
%ylabel('r^2'); plotPrefs;

% stats
statcnt = 1;
[pv,tbl,stats] = kruskalwallis(corr_all(:,mdlI),[],'off');
multc = multcompare(stats,'display','off');
stat(statcnt).type = 'kw test for model performance (all times)';
stat(statcnt).test = 'kruskalwallis';
stat(statcnt).p = pv;
stat(statcnt).median = median(corr_all(:,mdlI));
stat(statcnt).quartiles = prctile(corr_all(:,mdlI),[25 75]);
stat(statcnt).iqr = iqr(corr_all(:,mdlI));
stat(statcnt).n = size(corr_all(:,mdlI));
stat(statcnt).tbl = tbl;
stat(statcnt).stats = stats;
stat(statcnt).effectSizeMethod = 'eta2';
stat(statcnt).effectSize = (tbl{2,end-1}-size(corr_all(:,mdlI),2)+1) ...
    ./ (size(corr_all(:,mdlI),1) - size(corr_all(:,mdlI),2));
stat(statcnt).multc = multc;

% stats
statcnt = statcnt+1;
[pv,tbl,stats] = kruskalwallis(corr_not(:,mdlI),[],'off');
multc = multcompare(stats,'display','off');
stat(statcnt).type = 'kw test for model performance (after switch)';
stat(statcnt).test = 'kruskalwallis';
stat(statcnt).p = pv;
stat(statcnt).median = median(corr_not(:,mdlI));
stat(statcnt).quartiles = prctile(corr_not(:,mdlI),[25 75]);
stat(statcnt).iqr = iqr(corr_not(:,mdlI));
stat(statcnt).n = size(corr_not(:,mdlI));
stat(statcnt).tbl = tbl;
stat(statcnt).stats = stats;
stat(statcnt).effectSizeMethod = 'eta2';
stat(statcnt).effectSize = (tbl{2,end-1}-size(corr_not(:,mdlI),2)+1) ...
    ./ (size(corr_not(:,mdlI),1) - size(corr_not(:,mdlI),2));
stat(statcnt).multc = multc;


%% plot ln - glm correlations
ax3 = subplot(nrows,ncols,3); hold on
%for i = 1:3
%    x = gc(i+1,:); y = gc(1,:);
%    sh(i) = scatter(x,y,30,mdl_clite(2+i,:),'.');
    %end
[b,pv,xp,yp,ypci,mdl] = fitlmWrapper(gc(gc_metric,:)',gc(1,:)');
plotlmWrapper(xp,yp,ypci);
scatter(gc(gc_metric,:),gc(1,:),30,'k.');
axis tight;
plot([0 0],ylim,'k:');
plot(xlim,[0 0],'k:');
xlabel('GLM Gain Control');
ylabel('LN Gain Control');
title(sprintf('p = %g',pv(end))); plotPrefs;
axis square

% stats
statcnt = statcnt + 1;
stat(statcnt).type = sprintf('lm for gain control in %s vs lnmodel',...
                             ops.mdls{gc_metric+1});
stat(statcnt).test = 'lm';
stat(statcnt).p = pv;
stat(statcnt).n = size(gc(1,:));
stat(statcnt).tbl = mdl;
stat(statcnt).effectSize = mdl.Rsquared;
stat(statcnt).effectSizeMethod = 'r-squared';


%% plot gain control histogram
ax4 = subplot(nrows,ncols,4); hold on;
edges = linspace(min(gc(gc_metric,:)),max(gc(gc_metric,:)),20);
histogram(gc(gc_metric,:),edges,'facecolor','k');
mdn = median(gc(gc_metric,:));
plot([0 0],ylim,'k--','linewidth',1.5);
plot([mdn mdn],ylim,'r','linewidth',2);
[pv,~,stats] = signrank(gc(gc_metric,:));
xlabel('GLM Gain Control (high-low)'); ylabel('Cell Count');
title(sprintf('p = %g',pv)); plotPrefs;

% stats
statcnt = statcnt + 1;
[pv,~,stats] = signrank(gc(gc_metric,cell_include));
stat(statcnt).type = ['gain control different from 0, responsive ' ...
                    'neurons'];
stat(statcnt).test = 'signrank';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).median = median(gc(gc_metric,cell_include));
stat(statcnt).quartiles = prctile(gc(gc_metric,cell_include),[25 75]);
stat(statcnt).iqr = iqr(gc(gc_metric,cell_include));
stat(statcnt).n = sum(cell_include);
stat(statcnt).stats = stats;
stat(statcnt).effectSizeMethod = 'eta2';
stat(statcnt).effectSize = stats.zval / sqrt(sum(cell_include));

statcnt = statcnt + 1;
[pv,~,stats] = signrank(gc(gc_metric,:));
stat(statcnt).type = 'gain control different from 0, all neurons';
stat(statcnt).test = 'signrank';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).median = median(gc(gc_metric,:));
stat(statcnt).quartiles = prctile(gc(gc_metric,:),[25 75]);
stat(statcnt).iqr = iqr(gc(gc_metric,:));
stat(statcnt).n = length(gc);
stat(statcnt).stats = stats;
stat(statcnt).effectSizeMethod = 'eta2';
stat(statcnt).effectSize = stats.zval / sqrt(length(gc));


%%
% set neurons to include for adaptation analysis
%gc_cut = median(gc(gc_metric,:));
%include = gc(gc_metric,:) < gc_cut & ~any(true_tau==1,2)';
include = ~any(true_tau==1,2)' &  gc(gc_metric,:) < 0;

% plot average time course for included neurons
ax5 = subplot(nrows,ncols,5:6); hold on;
for i = 1:2
    
    y = mean(squeeze(glm_w(w_metric,ind+cind(i),include))');
    [params,mdl,tau] = fitExpGrid(xtimes(ind>=0),y(ind>=0),...
                                  [],[],[min(y)/2 -inf ops.period]);
    plot(xtimes,y,'.','color',c_lite(i,:));
    plot(xtimes(ind>=0),y(ind>=0),'.','color',c{i});
    plot(xf,mdl(params,xf),'color',c{i},'linewidth',1);
    
end
plot([0 0],ylim,'k--');
ylabel('w'); xlabel('Time (s)');
axis tight; plotPrefs;

ax6 = subplot(nrows,ncols,8:9); hold on;
for i = 1:2
    
    y = squeeze(glm_w(w_metric,ind+cind(i),include))';
    [params,mdl,tau] = fitExpGrid(xtimes(ind>=0),mean(y(:,ind>=0)),...
                                  [],[],[min(mean(y(:,ind>=0)))/2 -inf ops.period]);
    errorBars(xtimes,y,c_lite(i,:));
    %plot(xtimes(ind>=0),mean(y(:,ind>=0)),'-','color',c{i});
    plot(xf,mdl(params,xf),'color',c{i},'linewidth',1);
    
end
plot([0 0],ylim,'k--');
ylabel('w'); xlabel('Time (s)');
axis tight; plotPrefs;

% taus for included neurons
ax7 = subplot(nrows,ncols,7); hold on;
plotDist([1 2],true_tau(include,:),{'b','r'},true);
[pv,~,stats] = signrank(true_tau(include,1),true_tau(include,2));
title(sprintf('p = %g (n=%d)',pv,sum(include))); plotPrefs;
ylabel('\tau (s)');

% stats
statcnt = statcnt + 1;
stat(statcnt).type = ['gain tau between high and low contrast - ' ...
                    'neurons with GC < 0 and tau < 1']
stat(statcnt).test = 'signrank';
stat(statcnt).p = pv;
stat(statcnt).stats = stats;
stat(statcnt).median = median(true_tau(include,:));
stat(statcnt).quartiles = prctile(true_tau(include,:),[25 75]);
stat(statcnt).iqr = iqr(true_tau(include,:));
stat(statcnt).n = sum(include);
stat(statcnt).stats = stats;
stat(statcnt).effectSizeMethod = 'eta2';
stat(statcnt).effectSize = stats.zval / sqrt(sum(include));

ax4.Position([2 4]) = ax5.Position([2 4]);
ax7.Position([2 4]) = ax6.Position([2 4]);

saveFigPDF(f1,sz,sprintf('%s/_plots/_%sacute_2tau_summary.pdf',root,tag));


% plot example neuron
ex = 'K184_2021-04-21_13-48-12_052_181_mu';
mdl = 3;
% ex = 'K184_2021-04-21_13-48-12_028_076_mu';
exind = find(contains({r.cellID},ex));
file = dir(fullfile(resPath,sprintf('acute_%sdc%d*%s.mat',tag,dc(mdl),ex)));
load(fullfile(file.folder,file.name));

f2 = figure(2); set(f2,'Position',[0 0 900 900]); clf;
ops.mdls = models;
blks = 4; % number of blocks to repeat
twin = blks * ops.stimInfo.baseNoiseD;
ops.mdls = {'ln_static','ln-gc','glm-static','lohse','glm-symm','glm-asymm(1)','glm-asymm(2)'};
plot_acute_neuron_mdls(res,r(exind),u,s,ops,twin,[1 2 4 7],7);

saveFigPDF(f2,[900 900],sprintf('%s/_plots/_%sacute_2tau_example.pdf',root,tag));