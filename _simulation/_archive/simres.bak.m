clear all; close all;
addpath(genpath('~/chris-lab/code_general/'))

% run analysis for different parameters
gc = [0 .5 1];
correct = false;
simtype = 'sim3'
alpha = 0.95;

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
        mndenom(j,i) = min(abs(res.beta1+res.beta2));
        if isfield('res','gk')
            grange(j,i) = range(res.gk);
        else
            grange(j,i) = nan;
        end
        
        if correct
            if j == 1
                % load contrast and make design matrix
                cc = obj.mean_contrast./obj.contrast;
                cc(isinf(cc)) = 1;
                C = lagDesignMatrix(repmat(cc,1,length(ops.f)), ...
                                    length(ops.t));
            end

            
            %% eugenio fix
            beta1 = res.beta1;
            beta2 = res.beta2;
            gamma = beta2 .* C;
            dprod = (beta1 + beta2)' * (beta1+gamma);
            bbnorm = norm(beta1+beta2)^2;
            bgnorm = sum((beta1+gamma).^2,1);
            wnew = 2 * dprod ./ ...
                   (bbnorm - bgnorm + ...
                    sqrt((bbnorm+bgnorm).^2 - 4*(bbnorm*bgnorm- ...
                                                 dprod.^2)));
            wcorr = reshape(wnew',(ops.blockLength*2 + ops.pad)* ...
                            ops.fs,[])';
            wc(:,j,i) = mean(wcorr,1);
            errc(j,i) = immse(squeeze(wc(:,j,i)),squeeze(g(:,j,i)));
            
            

        end
                
        disp([i j]);
        
    end
    
end

f1 = figure(1); clf;
subplot(3,2,1)
cc = repmat([.66 .33 0],3,1)';
hold on;
for i = length(gc):-1:1;
    ph(i) = plot(mean(g(:,:,i),2),'--','Color',cc(i,:));
    errorbar(mean(w(:,:,i),2),sem(w(:,:,i)'),'-','Color',cc(i,:),'linewidth',1);
end
legend(ph,'g=1.0','g=0.5','g=0.0');
ylabel('w,g');
xlabel('Time (step)'); plotPrefs;
title(sprintf('alpha = %3.2f (mean)',alpha));
plotPrefs;

subplot(3,2,2)
hold on;
for i = length(gc):-1:1;
    ph(i) = plot(mean(g(:,:,i),2),'--','Color',cc(i,:));
    plot(median(w(:,:,i),2),'-','Color',cc(i,:),'linewidth',1);
end
legend(ph,'g=1.0','g=0.5','g=0.0');
ylabel('w,g');
xlabel('Time (step)'); plotPrefs;
title(sprintf('alpha = %3.2f (median)',alpha));
plotPrefs;

subplot(3,2,3)
plot(squeeze(w(:,:,3)),'k');
ylabel('w,g');
xlabel('Time (step)'); plotPrefs;

subplot(3,2,4)
ci = [ones(size(g,2),1) ones(size(g,2),1)*2 ones(size(g,2),1)*3];
cv = cc(ci(:),:);
scatter(log(err(:)),log(grange(:)),20,cv,'filled');
xlabel('log( mse )');
ylabel('logRange ( Gain Kernel )');
plotPrefs;

subplot(3,2,5);
D = log(pdist(w(:,:,3)','cosine'));
Z = squareform(D);
imagesc(Z); colorbar;
title('Dissimilarity of w');


saveFigPDF(f1,[800 900],sprintf('./_plots/_%s_res_alpha%03d.pdf', ...
                                simtype,alpha*100))

if correct
    f2 = figure(2); clf;
    subplot(3,2,1)
    cc = repmat([.66 .33 0],3,1)';
    hold on;
    for i = length(gc):-1:1;
        ph(i) = plot(mean(g(:,:,i),2),'--','Color',cc(i,:));
        errorbar(mean(wc(:,:,i),2),sem(wc(:,:,i)'),'-','Color',cc(i,:),'linewidth',1);
    end
    legend(ph,'g=1.0','g=0.5','g=0.0');
    ylabel('w,g');
    xlabel('Time (step)'); plotPrefs;
    title(sprintf('alpha = %3.2f (mean)',alpha));
    plotPrefs;

    subplot(3,2,2)
    hold on;
    for i = length(gc):-1:1;
        ph(i) = plot(mean(g(:,:,i),2),'--','Color',cc(i,:));
        plot(median(wc(:,:,i),2),'-','Color',cc(i,:),'linewidth',1);
    end
    legend(ph,'g=1.0','g=0.5','g=0.0');
    ylabel('w,g');
    xlabel('Time (step)'); plotPrefs;
    title(sprintf('alpha = %3.2f (median)',alpha));
    plotPrefs;

    subplot(3,2,3)
    plot(squeeze(wc(:,:,3)),'k');
    ylabel('w,g');
    xlabel('Time (step)'); plotPrefs;

    subplot(3,2,4)
    ci = [ones(size(g,2),1) ones(size(g,2),1)*2 ones(size(g,2),1)*3];
    cv = cc(ci(:),:);
    scatter(log(errc(:)),log(grange(:)),20,cv,'filled');
    xlabel('log( mse )');
    ylabel('logRange ( Gain Kernel )');
    plotPrefs;

    subplot(3,2,5);
    D = log(pdist(wc(:,:,3)','cosine'));
    Z = squareform(D);
    imagesc(Z); colorbar;
    title('Dissimilarity of w');

    saveFigPDF(f2,[800 900],sprintf('./_plots/_%s_res_gamma_alpha%03d.pdf', ...
                                    simtype,alpha*100))
end