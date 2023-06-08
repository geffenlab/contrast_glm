clear all; close all;
addpath(genpath('~/chris-lab/code_general/'))
addpath(genpath('../_functions/'));

%%
% this script analyzes glm fits to neural simulations, looking at
% effects of the degree of gain control (gc) and the amount of data
% used to fit (essentially, the number of unique noise exemplars, ex)

% run analysis for different parameters
gc = [1 .5 0 -.5 -1]
ex = [100 5];
nt = [5 100];
dc = [0 1 2];
tau= {'5-50'};
tag = '2tau-spline';
resDir = '~/data/gain_behavior/_simulations';

% load first file
f = dir(fullfile(resDir,sprintf('*%s*_gc%03d*ex%02d*nt%03d*tau%s*dc%d*.mat',tag,gc(1)*10,ex(1),nt(1),tau{1},dc(1))));
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
h = colorbar; colormap(ax1,flipud(cbrewer2('PiYG')));
caxis([clims(1) clims(end)]);
plotPrefs;
saveFigPDF(ax1,[400 400],'./_plots/_strf_example.pdf');


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
                                     sprintf('sim_%s*_gc%03d_ex%02d_nt%03d_tau%se-2*dc%d*.mat',tag,...
                                             gc(j)*10,ex(i),nt(i),tau{ii},dc(m))));
                
                if length(files) < 100
                    error('Not enough files');
                end
                
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
colors_lite = flipud([.75 0 .75; ...
                    .5 0 .5; ...
                    .5 .5 .5; ...
                    0 .5 .5; ...
                    0 .75 .75]+.25);
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

    for ii = 1:length(ex)
        
        exI = ii;

        for j = 1:length(tau)
            
            
            %% stats
            strf = squeeze(mean(squeeze(mdl(1).strf(:,i,exI,1,:,:)),1,'omitnan'));
            c = mean(reshape(neuron.ops.contrast,ns,[])',1);
            g = squeeze(mean(n.g(:,i,exI,j,:),1));
            y = squeeze(mean(n.y(:,i,exI,j,:),1));
            x = 1:length(y);
            ld = reshape(squeeze(n.ld(exI,i,j,:)),ns,[])';
            for k = 1:numel(mdl)
                p_full(k,:) = squeeze(mean(mdl(k).pred(:,i,exI,j,:),1));
                w_full(k,:) = squeeze(mean(mdl(k).w(:,i,exI,j,:),1));
                p_std(k,:) = squeeze(std(mdl(k).pred(:,i,exI,j,:),1));
                w_std(k,:) = squeeze(std(mdl(k).w(:,i,exI,j,:),1));
                w_all(k,:,:) = squeeze(mdl(k).w(:,i,exI,j,:));
            end
            
            
            %% plot overlaid figures
            figure(f22);
            subplot(2,2,1+(ii-1)*2); hold on;
            plot(log(y),'color',colors(i,:),'linewidth',lw(j),'linestyle','-');
            plot(log(p_full(3,:)),'color',[0 0 0 .5],'linestyle','--','linewidth',1);
            ylabel('log(spike count)'); plotPrefs;
            
            subplot(2,2,2+(ii-1)*2); hold on;
            %errorbar(w_full(3,:),w_std(3,:),'color',[colors(i,:) ...
%                    .25],'linestyle','--','linewidth',.5);
            [ph,lh] = patchErrorBars(1:length(w_full),squeeze(w_all(3,:,:)),colors(i,:),'prctile');
            lh.LineStyle = '--';
            lh.Color(4) = .5;
            ph.FaceAlpha = .2;
            plot(g,'color',colors(i,:),'linewidth',1.5,'linestyle','-');
            ylabel('Gain'); plotPrefs;
            
            
        end
        
    end
    
    %fn = sprintf('./_plots/_%s_3mdls_gc%d.pdf',tag,gc(i)*10);
    %saveFigPDF(fh(i),sz,fn);
    
end

saveFigPDF(f22,[600 400],'./_plots/_ex_comp_2tau.pdf');









