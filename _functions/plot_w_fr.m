function [contrast dglm w_mat psth] = plot_w_fr(d,pdat,include,ops)

% fix up formatting of cell ids
glmids = vertcat(d.cellID{include});
glmids = cellstr(glmids(:,16:end-3));
glmids = regexprep(glmids,'-','');
pids = vertcat(pdat.res.cellID{:});
pids = cellstr(pids(:,[1:17 19:end-3]));
pids = regexprep(pids,'_','');

% match
[lia,lib] = ismember(glmids,pids);

% extract
w_mat = d.w_new(include,:);
w_mat = w_mat(lia,:);
psth = pdat.res.mPSTH(lib(lib>0),:);
psthz = normalize(psth,2);
psthv = pdat.res.vPSTH(lib(lib>0),:);
contrast = d.contrastI(include);
contrast = contrast(lia);
wt = (0:size(w_mat,2)-1) / ops.fs;
pt = pdat.ops.time;
dglm = mean(w_mat(:,140:160),2) - mean(w_mat(:,100:120),2);
dglm(contrast==1) = -dglm(contrast==1);

nrows = 8;
ncols = 5;
cols = {'b','r'};
cols_lite = {[.5 .5 1]; [1 .5 .5]};


%figure;
%[~,sortI] = sortrows([contrast',dglm]);
%imagesc(w_mat(sortI,:),[.5 1.5]);


for i = 1:2
    
    cI = contrast == (i-1);
    
    % w and dw for this contrast
    w = w_mat(cI,:);
    dw = dglm(cI);
    
    % quantiles
    prct = [0 .05 .25 .75 .95 1];
    q = quantile(dw,prct);
    edges = linspace(-1,.5,30);
    
    % for each quantile, make a plot
    qI = zeros(size(dw));
    for j = 1:length(q)-1
        qI(dw>q(j)&dw<q(j+1)) = j;
        
        subplot(nrows,ncols,(4*ncols)*(i-1)+j); hold on;
        histogram(dw,edges,'FaceColor',cols_lite{i},'FaceAlpha',.1, ...
                  'LineStyle','--');
        histogram(dw(qI==j),edges,'FaceColor',cols_lite{i});
        plot([0 0],ylim,'k','LineWidth',1);
        plot([q(j) q(j)],ylim,'Color',cols{i});
        plot([q(j+1) q(j+1)],ylim,'Color',cols{i});
        plotPrefs; 
        if i == 1; 
            title(sprintf('%d-%d percentile\nw(lo) - w(hi)', ...
                          prct(j)*100,prct(j+1)*100));
        end
        
        subplot(nrows,ncols,(4.*ncols*(i-1))+[j+5 j+10]); hold on;
        imagesc(wt,1:sum(qI==j),w(qI==j,:),[.5 1.5]); axis tight;
        xlim([2 4]); plotPrefs;
        
        subplot(nrows,ncols,(4.*ncols*(i-1))+j+15);
        plotPrefs;
        %yyaxis right;
        %plot(pt,SmoothGaus(median(psthv(qI==j,:)),2),...
%     'Color',[.5 .5 .5]);
%if j == 1; ylabel('var(FR)'); end
%        xlim([2 4]);
%        yyaxis left;
        plot(pt,SmoothGaus(median(psthz(qI==j,:)),1),'color', ...
             cols{i},'linewidth',1);
        xlim([2 4]); ylim([-.5 1]);
        if j == 1; ylabel('FR (z)'); end
        
        
    end
    
end