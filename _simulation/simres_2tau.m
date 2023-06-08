clear all; close all;
addpath(genpath('~/chris-lab/code_general/'))
addpath(genpath('../_functions/'));

%%
% this script analyzes glm fits to neural simulations, looking at
% effects of the degree of gain control (gc) and the amount of data
% used to fit (essentially, the number of unique noise exemplars, ex)

% run analysis for different parameters
gc = [1 .5 0 -.5 -1]
ex = [100];
nt = [5];
dc = [0 1 2];
tau= {'5-50','50-50','50-5'};
tag = '2tau-spline';
resDir = '~/data/gain_behavior/_simulations';

% load first file
f = dir(fullfile(resDir,sprintf('*%s_gc%03d*ex%02d*nt%03d*tau%s*dc%d*.mat',tag,gc(1)*10,ex(1),nt(1),tau{1},dc(1))));
load(fullfile(f(1).folder,f(1).name))

% preallocate
ntf = length(ops.t); nff = length(ops.f);
ns = (ops.blockLength*2 + ops.pad)*ops.fs;

n.y = nan(length(f),length(gc),length(ex),length(tau),ns);
n.g = n.y;

strf = neuron.ops.beta;
ax1 = figure;
clims = linspace(-max(abs(strf(:))),max(abs(strf(:))),1000);
cmap = zeroCMap(clims,0);
imagesc(ops.t,ops.f/1000,strf);
h = colorbar; colormap(ax1,cmap);
caxis([clims(1) clims(end)]);
plotPrefs;

for i = 1:3
    
    f = dir(fullfile(resDir,sprintf('*%s_gc%03d*ex%02d*nt%03d*tau%se-2*dc%d*.mat',tag,gc(1)*10,ex(1),nt(1),tau{1},dc(i))));
    load(fullfile(f(1).folder,f(1).name))
    
    mdl(i).w = n.y;
    mdl(i).b0 = nan(length(f),length(gc),length(ex),length(tau));
    mdl(i).b1 = mdl(i).b0;
    mdl(i).b2 = nan(length(f),length(gc),length(ex),length(tau),length(res.beta2));
    mdl(i).b3 = nan(length(f),length(gc),length(ex),length(tau),length(res.beta3));
    mdl(i).pred = n.y;
    mdl(i).strf = nan(length(f),length(gc),length(ex),length(tau),nff,ntf);

    if isfield(res,'basis')
        if dc(i) == 0
            % symmetric kernels
            mdl(i).b2_kernel = nan(length(f),length(gc),length(ex),length(tau),ops.gain_lags);
            mdl(i).b3_kernel = nan(length(f),length(gc),length(ex),length(tau),ops.gain_lags);
        elseif dc(i) == 1
            % split gain kernel, symmetric contrast kernel
            mdl(i).b2_kernel = nan(length(f),length(gc),length(ex),length(tau),ops.gain_lags,2);
            mdl(i).b3_kernel = nan(length(f),length(gc),length(ex),length(tau),ops.gain_lags);
        elseif dc(i) == 2
            % split gain kernel and contrast kernel
            mdl(i).b2_kernel = nan(length(f),length(gc),length(ex),length(tau),ops.gain_lags,2);
            mdl(i).b3_kernel = nan(length(f),length(gc),length(ex),length(tau),ops.gain_lags,2);
        end
    end
    
end


% for each exemplar type
for i = 1:length(ex)
    
    % for each gain control value
    for j = 1:length(gc)
        
        % for each transition
        for ii = 1:length(tau)
            
            % for each model
            for m = 1:length(dc)
                
                
                % build file list for this parameter set
                files = dir(fullfile(resDir,...
                                     sprintf('sim_%s_gc%03d_ex%02d_nt%03d_tau%se-2*dc%d*.mat',tag,...
                                             gc(j)*10,ex(i),nt(i),tau{ii},dc(m))));
                
                for k = 1:length(files)
                    
                    load(fullfile(files(k).folder,files(k).name), ...
                         'neuron','res','ops');
                    
                    %% neuron stuff
                    n.g(k,j,i,ii,:) = neuron.g(1:ns);
                    n.y(k,j,i,ii,:) = mean(reshape(neuron.y,ns,[])',1);
                    
                    
                    %% full model params
                    mdl(m).b0(k,j,i,ii) = res.beta0;
                    mdl(m).b1(k,j,i,ii) = res.beta1;
                    mdl(m).b2(k,j,i,ii,:) = res.beta2;
                    mdl(m).b3(k,j,i,ii,:) = res.beta3;
                    mdl(m).w(k,j,i,ii,:) = mean(reshape(res.w,ns,[])',1);
                    mdl(m).pred(k,j,i,ii,:) = mean(reshape(res.gain_fit.pred,ns,[])',1);
                    mdl(m).strf(k,j,i,ii,:,:) = res.strf;
                    
                    if isfield(res,'basis')
                        if dc(m) == 0
                            mdl(m).b2_kernel(k,j,i,ii,:) = res.basis*res.beta2;
                            mdl(m).b3_kernel(k,j,i,ii,:) = res.basis*res.beta3;
                        elseif dc(m) == 1
                            mdl(m).b2_kernel(k,j,i,ii,:,1) = res.basis*res.beta2(1:length(res.beta2)/2);
                            mdl(m).b2_kernel(k,j,i,ii,:,2) = res.basis*res.beta2(length(res.beta2)/2+1:end);
                            mdl(m).b3_kernel(k,j,i,ii,:) = res.basis*res.beta3;
                        elseif dc(m) == 2
                            mdl(m).b2_kernel(k,j,i,ii,:,1) = res.basis*res.beta2(1:length(res.beta2)/2);
                            mdl(m).b2_kernel(k,j,i,ii,:,2) = res.basis*res.beta2(length(res.beta2)/2+1:end);
                            mdl(m).b3_kernel(k,j,i,ii,:,1) = res.basis*res.beta3(1:length(res.beta3)/2);
                            mdl(m).b3_kernel(k,j,i,ii,:,2) = res.basis*res.beta3(length(res.beta3)/2+1:end);
                        end
                    end
                    

                    
                    
                    fprintf('%s %d exemplars, %03.1f gain control, tau %s, dc %d, run %03d\n',...
                            tag,ex(i),gc(j),tau{ii},dc(m),k);
                    
                end
                
                n.ld(i,j,ii,:) = neuron.h;
                
            end
                
        end
                
    end
    
end

keyboard


colors = flipud([.75 0 .75; ...
                 .5 0 .5; ...
                 .5 .5 .5; ...
                 0 .5 .5; ...
                 0 .75 .75]);
ls = {':','--','-.'};
ms = {'.','x','^','*'};
lw = [1 1 1];
nrows = 5;
ncols = 14;
my = mean(log(n.y),1);
yplim = [floor(min(my(~isinf(my(:))))) ceil(max(my(:)))];
cc = [0 0 1; 1 0 0];
cc_lite = [.5 .5 1; 1 .5 .5];
exI = 1;
sz = [2400 400];
nrows = 3; ncols = 11;
mdls = {'Symm. Model','Asymm. Model (1)','Asymm. Model (2)'};
mdl_col = [27,158,119;...
           217,95,2;...
           117,112,179] ./ 256;
mdl_col_lite = [102,194,165;...
                252,141,98;...
                141,160,203] ./ 256;

% figure with overlaid
f22 = figure(22); clf;

for i = 1:length(gc)
    fh(i) = figure(i); clf;
    set(fh(i),'Position',[0 0 sz]);

    for j = 1:length(tau)
        
        figure(fh(i));
        
        %% stats
        strf = squeeze(mean(squeeze(mdl(1).strf(:,i,exI,:,:,:)),[1 2]));
        c = mean(reshape(neuron.ops.contrast,ns,[])',1);
        g = squeeze(mean(n.g(:,i,exI,j,:),1));
        y = squeeze(mean(n.y(:,i,exI,j,:),1));
        x = 1:length(y);
        ld = reshape(squeeze(n.ld(exI,i,j,:)),ns,[])';
        for k = 1:numel(mdl)
            p_full(k,:) = squeeze(mean(mdl(k).pred(:,i,exI,j,:),1));
            w_full(k,:) = squeeze(mean(mdl(k).w(:,i,exI,j,:),1));
        end
        
        
        %% strf
        ax1 = subplot(nrows,ncols,1+(j-1)*ncols);
        clims = linspace(-max(abs(strf(:))),max(abs(strf(:))),1000);
        cmap = zeroCMap(clims,0);
        imagesc(ops.t,ops.f/1000,strf);
        h = colorbar; colormap(ax1,cmap);
        caxis([clims(1) clims(end)]);
        plotPrefs;
        if j == 1; title('STRF (dc=0)'); end;

        
        %% linear drive
        ax2 = subplot(nrows,ncols,2+(j-1)*ncols);
        clims = linspace(-max(abs(ld(:))),max(abs(ld(:))),1000);
        imagesc(ld); h = colorbar; colormap(ax2,cmap);
        caxis([clims(1) clims(end)]);
        plotPrefs;

        
        %% firing rate/predictions
        subplot(nrows,ncols,[3 4]+(j-1)*ncols); hold on; clear pp;
        %pp(4) = plot(c,'-','color','k','linewidth',.5);
        %plot(x(c==min(c)),c(c==min(c)),'.','color',cc_lite(1,:));
        %plot(x(c==max(c)),c(c==max(c)),'.','color',cc_lite(2,:));
        pp(1) = plot(y,'color',colors(i,:),'linewidth',lw(j),'linestyle','-');
        for k = 1:length(mdl)
            pp(k+1) = plot(p_full(k,:),'color','k','linestyle',ls{k},'linewidth',1);
        end
        YL = ylim;
        plot([1 80],[YL(2) YL(2)],'-','color',cc_lite(1,:),'linewidth',2);
        plot([81 160],[YL(2) YL(2)],'-','color',cc_lite(2,:),'linewidth',2);
        ylabel('spike count'); plotPrefs;
        legend(pp,['Spikes',mdls],'location','northeastoutside');

        
        %% gain/w
        subplot(nrows,ncols,[5 6]+(j-1)*ncols); hold on; clear pp;
        plot([1 80],[1.6 1.6],'-','color',cc_lite(1,:),'linewidth',2);
        plot([81 160],[1.6 1.6],'-','color',cc_lite(2,:),'linewidth',2);
        pp(1) = plot(g,'color',colors(i,:),'linewidth',lw(j));
        for k = 1:length(mdl)
            pp(k+1) = plot(w_full(k,:),'color','k','linestyle',ls{k},'linewidth',1);
        end        
        legend(pp,['True Gain',mdls],'location','northeastoutside');
        ylabel('Gain'); plotPrefs;
        
        
        %% coefficients
        ax = subplot(nrows,ncols,[7 8]+(j-1)*ncols); hold on;
        offset = 0;
        x1 = [0 10 20 30+2*size(mdl(1).b2,5)];
        for k = 3:-1:1
            clear xx;
            xx{1} = x1(1);
            xx{2} = x1(2);
            xx{3} = x1(3):x1(3)+size(mdl(k).b2,5)-1;
            xx{4} = x1(4):x1(4)+size(mdl(k).b3,5)-1;
            for ii = 1:length(xx)
                lbl = sprintf('b%d',ii-1);
                %ee(1) = errorbar(xx{ii},mean(squeeze(mdl(k).(lbl)(:,i,exI,j,:)),1),...
%                 std(squeeze(mdl(k).(lbl)(:,i,exI,j,:)),1),...
%                                'color',mdl_col_lite(k,:));
                ee(k) = plot(xx{ii}+offset*(k-1),mean(squeeze(mdl(k).(lbl)(:,i,exI,j,:)),1),...
                             'color',mdl_col(k,:),'marker','.');
            end
            
        end
        plotPrefs; symlog(ax,'y',-1.5); set(ax,['y','MinorGrid'],'off');
        legend(ee,mdls,'location','northeastoutside');
        xlim([-5 x1(4)+size(mdl(3).b3,5)+5+offset*length(mdls)-1]);
        set(gca,'xtick',x1,'xticklabels',{'b0','b1','b2','b3'});
        ylabel('Coefficient');
        
        
        %% kernels for gain and contrast if splines are used
        for k = 1:length(mdl)
            if isfield(mdl(k),'b2_kernel')
                subplot(nrows,ncols,9+(j-1)*ncols+(k-1)); hold on
                t = (0:ops.gain_lags-1)/ops.fs;
                if dc(k) == 0
                    plot(t,mean(squeeze(mdl(k).b3_kernel(:,i,exI,j,:))),'k--','linewidth',1);
                    plot(t,mean(squeeze(mdl(k).b2_kernel(:,i,exI,j,:))),'k','linewidth',1);
                    legend('Contrast Kernel','Gain Kernel');
                elseif dc(k) == 1
                    plot(t,mean(squeeze(mdl(k).b3_kernel(:,i,exI,j,:))),'k--','linewidth',1);
                    plot(t,mean(squeeze(mdl(k).b2_kernel(:,i,exI,j,:,1))),'b','linewidth',1);
                    plot(t,mean(squeeze(mdl(k).b2_kernel(:,i,exI,j,:,2))),'r','linewidth',1);
                elseif dc(k) == 2
                    plot(t,mean(squeeze(mdl(k).b3_kernel(:,i,exI,j,:,1))),'b--','linewidth',1);
                    plot(t,mean(squeeze(mdl(k).b3_kernel(:,i,exI,j,:,2))),'r--','linewidth',1);
                    plot(t,mean(squeeze(mdl(k).b2_kernel(:,i,exI,j,:,1))),'b','linewidth',1);
                    plot(t,mean(squeeze(mdl(k).b2_kernel(:,i,exI,j,:,2))),'r','linewidth',1);
                end
            end
            xlabel('Time (s)'); ylabel('Kernel Weight');
            plotPrefs; title(['Kernels (' mdls{k} ')']);
        end
        
        
        %% plot overlaid figures
        figure(f22);
        subplot(3,2,1+(j-1)*2); hold on;
        plot(log(y),'color',colors(i,:),'linewidth',lw(j),'linestyle','-');
        plot(log(p_full(3,:)),'color',[0 0 0 .5],'linestyle',ls{3},'linewidth',1);
        ylabel('log(spike count)'); plotPrefs;
        
        subplot(3,2,2+(j-1)*2); hold on;
        plot(g,'color',colors(i,:),'linewidth',lw(j),'linestyle','-');
        plot(w_full(3,:),'color',[0 0 0 .5],'linestyle',ls{3},'linewidth',1);
        ylabel('Gain'); plotPrefs;
        
                
    end
    
    fn = sprintf('./_plots/_%s_3mdls_gc%d.pdf',tag,gc(i)*10);
    saveFigPDF(fh(i),sz,fn);
    
end

saveFigPDF(f22,[600 400],'./_plots/_gc_comp_2tau.pdf');









