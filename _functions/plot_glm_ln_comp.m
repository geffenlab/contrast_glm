function plot_glm_ln_comp(d,g,lut,include,ops)
        
% get gain at steady states from contrast glm
glm_g = [mean(d.w_new(include,100:120),2) ...
            mean(d.w_new(include,140:160),2)];


% match cells by firing rate (not ideal but seems to work)
a = g.t.mfr;
b = lut.FR(include);
c = d.contrastI(include);
[lib,lia] = ismember(b,a);

% ln values
ln_gain = g.t.glm_ahat(lia(lia>0),:,3);
ln_contrast = g.t.contrastI(lia(lia>0));

% glm values
glm_gain = glm_g(lib,:);
glm_contrast = c(lib)';

% joint contrast
contrast = glm_contrast .* ln_contrast;

% gain changes (multiply by -contrast to get proper direction)
dglm = glm_gain(:,2) - glm_gain(:,1);
dglm(contrast==1) = -dglm(contrast==1);
dln = (ln_gain(:,2) - ln_gain(:,1));
dln(contrast==1) = -dln(contrast==1);


bins = 50;
edges = linspace(0,.6,bins);
subplot(2,3,1); hold on;
histogram(ln_gain(contrast==1,1),edges,'FaceColor','b');
histogram(ln_gain(contrast==1,2),edges,'FaceColor','r');
histogram(ln_gain(contrast==0,1),edges,'FaceColor','r','FaceAlpha',.1,'LineStyle','--');
histogram(ln_gain(contrast==0,2),edges,'FaceColor','b','FaceAlpha',.1,'LineStyle','--');
legend('L->H (Adaptation Period)','L->H (Target Period)',...
       'H->L (Adaptation Period)','H->L (Target Period)');
xlabel('LN Model Gain'); ylabel('Cell Count'); plotPrefs;

subplot(2,3,2); hold on;
edges = linspace(.4,1.6,bins);
histogram(glm_gain(contrast==1,1),edges,'FaceColor','b');
histogram(glm_gain(contrast==1,2),edges,'FaceColor','r');
histogram(glm_gain(contrast==0,1),edges,'FaceColor','r','FaceAlpha',.1,'LineStyle','--');
histogram(glm_gain(contrast==0,2),edges,'FaceColor','b','FaceAlpha',.1,'LineStyle','--');
xlabel('GLM Gain'); ylabel('Cell Count'); plotPrefs;

subplot(2,3,3); hold on;
plot([0 0],[-1 1],'k--');
plot([-1 1],[0 0],'k--');
sh(2) = scatter(dln(glm_contrast==0),dglm(contrast==0),20,...
                [.5 .5 .5],'.');
sh(1) = scatter(dln(glm_contrast==1),dglm(contrast==1),20,'k.');
xlim([-1 1]); ylim([-1 1]); plotPrefs;
xlabel('LN gain diff (lo-hi)');
ylabel('GLM gain diff (lo-hi)');
legend(sh,'L->H','H->L');

subplot(2,3,4); hold on;
edges = linspace(-1,1,bins);
hh(1) = histogram(dln(contrast==1),edges,'facecolor','k');
hh(2) = histogram(dln(contrast==0),edges,'facecolor','k','FaceAlpha',.1, ...
          'LineStyle','--');
plot([0 0],ylim,'k','linewidth',1);
xlabel('LN \Deltagain (lo-hi)'); ylabel('Cell Count'); plotPrefs;
legend(hh,'L->H','H->L');

subplot(2,3,5); hold on;
histogram(dglm(contrast==1),edges,'facecolor','k')
histogram(dglm(contrast==0),edges,'facecolor','k','FaceAlpha',.1, ...
          'LineStyle','--')
plot([0 0],ylim,'k','linewidth',1);
xlabel('GLM \Deltagain (lo-hi)'); ylabel('Cell Count'); plotPrefs;

% visualize all good fits sorted by contrast (look only around the
% transition)
subplot(2,3,6);hold on;
[~,sortI] = sortrows([d.contrastI(include)' glm_g(:,1)-glm_g(:,2)]);
w_mat = d.w_new(include,100:end);
tv = (100:size(d.w_new,2)) ./ ops.fs;
imagesc(tv,1:length(sortI),w_mat(sortI,:),[.5 1.5]);
plot(xlim,[sum(d.contrastI(include)==0) sum(d.contrastI(include)==0)], ...
     'r','linewidth',1);
set(gca,'YTick',[600 1600]);
set(gca,'YTickLabels',{'H->L','L->H'})
colorbar; plotPrefs; axis tight;
xlabel('Time (s)'); ylabel('Cells (sort by contrast, \Deltagain)')
title('GLM Gain Estimates (w)')
