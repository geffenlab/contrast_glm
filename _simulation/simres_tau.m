clear all; close all;
addpath(genpath('~/chris-lab/code_general/'))
addpath(genpath('../_functions/'));

%%
% this script analyzes glm fits to neural simulations, looking at
% effects of the degree of gain control (gc) and the amount of data
% used to fit (essentially, the number of unique noise exemplars, ex)

% run analysis for different parameters
gc = [-1 -.5 0 .5 1]
ex = [5 100]
nt = [100 5]
tau= [10 .5 .2 .1]
tag = 'lohi-notrans_';
resDir = fullfile('_res');

% load first file
f = dir(fullfile(resDir,sprintf('sim_%sgc%03d*ex%02d*nt%03d*tau%03d*.mat',tag,gc(1)*10,ex(1),nt(1),tau(1)*10)));
load(fullfile(f(1).folder,f(1).name))

% preallocate
ntf = length(ops.t); nff = length(ops.f);
ns = (ops.blockLength*2 + ops.pad)*ops.fs;
w = nan(length(f),length(gc),length(ex),length(tau),ns);
g = w;
strf = nan(length(f),length(gc),length(ex),length(tau),length(res.strf));
b0 = nan(length(f),length(gc),length(ex),length(tau));
b1 = b0;
b2 = nan(length(f),length(gc),length(ex),length(tau),length(res.beta2));
b3 = b2;
y = nan(length(f),length(gc),length(ex),length(tau),ns);
ld = nan(length(gc),length(ex),ns*ops.ntrials*ops.exemplars);

% for each exemplar type
for i = 1:length(ex)
    
    % for each gain control value
    for j = 1:length(gc)
        
        for ii = 1:length(tau)
            
                    
            % build file list for this parameter set
            files = dir(fullfile(resDir,...
                                 sprintf('sim_%sgc%03d*ex%02d*nt%03d*tau%03d*.mat',tag,...
                                         gc(j)*10,ex(i),nt(i),tau(ii)*10)));
            
            for k = 1:length(files)
                
                load(fullfile(files(k).folder,files(k).name),'neuron','res','ops');
                
                % compile time course, true gain, compute error, range etc
                g(k,j,i,ii,:) = neuron.g(1:ns);
                b0(k,j,i,ii) = res.beta0;
                b1(k,j,i,ii) = res.beta1;
                b2(k,j,i,ii,:) = res.beta2;
                b3(k,j,i,ii,:) = res.beta3;
                strf(k,j,i,ii,:) = res.strf;
                y(k,j,i,ii,:) = mean(reshape(neuron.y,ns,[])',1);
                
                % gain estimate (w)
                w(k,j,i,ii,:) = mean(res.wm,1);
                
                fprintf('%s %d exemplars, %03.1f gain control, tau %03.1f, run %03d\n',...
                        tag,ex(i),gc(j),tau(ii),k);
                
            end
            
        end
        
        ld(j,i,:) = neuron.h;
        
    end
    

end



for k = 1:length(ex)

    % setup plot
    f1 = figure(1); clf;
    colors = [.75 0 .75; ...
              .5 0 .5; ...
              .5 .5 .5; ...
              0 .5 .5; ...
              0 .75 .75];
    ls = {'-','--',':','-.'};
    ms = {'.','x','^','*'};
    lw = [.25 .75 1 1.25];
    for i = 1:length(gc)
        lg{i} = sprintf('g=%3.1f',gc(i));
    end
    nrows = 5;
    ncols = 14;
    exI = k;
    my = mean(log(y),1);
    yplim = [floor(min(my(~isinf(my(:))))) ceil(max(my(:)))];

    
    % plot each gain condition over rows
    for i = 1:length(gc)
        
        subplot(nrows,ncols,(ncols)*(i-1)+[1]); hold on;
        mstrf = squeeze(mean(squeeze(strf(:,i,exI,:,2:end)),[1 2]));
        imagesc(ops.t,ops.f/1000,reshape(mstrf,nff,[])); colorbar;
        axis tight; plotPrefs;
        if i == 1; title('STRF Estimate'); end
        if i == length(gc); xlabel('Lag (s)'); ylabel('Freq. (kHz)'); end

        
        subplot(nrows,ncols,(ncols)*(i-1)+[2]); hold on;
        imagesc(reshape(ld(i,exI,:),ns,[])'); colorbar;
        if i == 1; title(sprintf('%d FROZEN NOISE PATTERNS\n\nLinear Drive',ex(exI))); end
        if j == length(gc); xlabel('Time (step)'), ylabel('Trial'); end
        plotPrefs; axis tight;
        
        for j = 1:length(tau)
            
            % y
            subplot(nrows,ncols,(ncols)*(i-1)+[3:4]); hold on;
            plot(1:ns,mean(log(squeeze(y(:,i,exI,j,:)))),'color',colors(i,:),'LineWidth',lw(j));
            plotPrefs; axis tight;
            ylim(yplim);
            if i == 1; title('Spike Rate'); end;
            if i == length(gc); xlabel('Time (step)'); ylabel('log(y)');  end
            
            % y (zoom)
            subplot(nrows,ncols,(ncols)*(i-1)+[5]); hold on;
            lnh(j) = plot(1:ns,mean(log(squeeze(y(:,i,exI,j,:)))),'color',colors(i,:),'LineWidth',lw(j));
            plotPrefs; axis tight;
            ylim(yplim); xlim([75 120]);
            if i == 1; title('Spike Rate (zoom)'); end
            if i == length(gc); xlabel('Time (step)'); end
            
            % w
            subplot(nrows,ncols,(ncols)*(i-1)+[6 7]); hold on;
            plot(1:ns,mean(squeeze(w(:,i,exI,j,:))),'color',colors(i,:),'LineWidth',lw(j));
            plot(squeeze(mean(g(:,i,exI,j,:),1)),'--','Color',colors(i,:),'LineWidth',lw(j));
            plotPrefs; axis tight;
            ylim([.5 1.5]);
            if i == 1; title('Gain'); end;
            if i == length(gc); xlabel('Time (step)'); ylabel('w'); end

            
            % w (zoom)
            subplot(nrows,ncols,(ncols)*(i-1)+[8]); hold on;
            plot(1:ns,mean(squeeze(w(:,i,exI,j,:))),'color',colors(i,:),'LineWidth',lw(j));
            plot(squeeze(mean(g(:,i,exI,j,:),1)),'--','Color',colors(i,:),'LineWidth',lw(j));
            plotPrefs; axis tight;
            ylim([.5 1.5]); xlim([75 120]);
            if i == 1; title('Gain Zoom'); end;
            if i == length(gc); xlabel('Time (step)'); end
            
            % betas
            subplot(nrows,ncols,(ncols)*(i-1)+[9]); hold on;
            errorbar(j,mean(b0(:,i,exI,j)),sem(b0(:,i,exI,j)),...
                     'color',colors(i,:),'linewidth',lw(j),'marker','o');
            plotPrefs; axis tight;
            set(gca,'xtick',1:length(tau));
            set(gca,'xticklabels',num2str(tau'));
            if i == 1; title('Baseline Rate'); end;
            if i == length(gc); xlabel('\tau'); ylabel('\beta_0'); end
            xlim([0 5]);
            
            subplot(nrows,ncols,(ncols)*(i-1)+[10]); hold on;
            errorbar(j,mean(b1(:,i,exI,j)),sem(b1(:,i,exI,j)),...
                     'color',colors(i,:),'linewidth',lw(j),'marker','o');
            plotPrefs; axis tight;
            set(gca,'xtick',1:length(tau));
            set(gca,'xticklabels',num2str(tau'));
            if i == 1; title('STRF Predictor'); end;
            if i == length(gc); xlabel('\tau'); ylabel('\beta_1'); end
            xlim([0 5]);
            
            subplot(nrows,ncols,(ncols)*(i-1)+[11:12]); hold on;
            bb2 = squeeze(b2(:,i,exI,j,:));
            [pph plh] = patchErrorBars(5:5+size(bb2,2)-1,bb2,colors(i,:));
            plh.LineWidth = lw(j);
            pph.FaceAlpha = .1;
            plotPrefs; axis tight;
            if i == 1; title('Gain Predictors'); end;
            if i == length(gc); xlabel('Lag (step)'); ylabel('\beta_2'); end
            
            subplot(nrows,ncols,(ncols)*(i-1)+[13:14]); hold on;
            bb3 = squeeze(b3(:,i,exI,j,:));
            [pph plh] = patchErrorBars(5+size(bb2,2)+1:5+2*size(bb3,2),bb3,colors(i,:));
            plh.LineWidth = lw(j);
            pph.FaceAlpha = .1;
            plotPrefs; axis tight;
            if i == 1; title('Contrast Predictor'); end;
            if i == length(gc); xlabel('Lag (step)'); ylabel('\beta_3'); end

        end
        
    end
    

    saveFigPDF(f1,[2400 600],sprintf(['./_plots/' ...
                        '_res_2step_%s%s_ex%03d.pdf'],tag,'tau',ex(k)));

end
