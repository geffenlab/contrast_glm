function [res,ops] = fitLN(res,y0,cc0,ops,X)

if isstruct(res) & isfield(res,'strf_fit_cv')
    
    if ~exist('X','var') | isempty(X)
        % if res is a struct with stf_fit_cv, use the previous
        % crossvalidation runs for the strf prediction
        
        % reconstruct prediction
        x = nan(size(ops.include));
        for i = 1:ops.cvfolds
            x(ops.cv_sample ~= i) = res.strf_fit_cv(i).pred;
        end
    end
    
elseif ~isstruct(res) & any(size(res)==1)
    % if res is a vector, set up data for cross validation
    ops.cv_trials = crossvalind('kfold',ops.order_r(:,1),ops.cvfolds);
    ops.cv_sample = cleanResample(ops.cv_trials',ops.blockLength,ops.period)';
    x = res;
    
end

% create a mask for transition times
tmp = [ones(1,ops.fs) zeros(1,ops.fs*(ops.blockLength-1))];
trans = repmat(tmp,1,length(ops.include)/length(tmp))';

% contrast values for indexing
ch = unique(cc0,'stable');

% cross-validated fits for each set of predictions
fprintf(sprintf('Fitting LN models with %d folds: ',ops.cvfolds)); tic;
for i = 1:ops.cvfolds
    
    % include sample trials, exclude transition
    I = ops.cv_sample ~= i & ~trans;
        
    % get contrast index for each sample
    cc = cc0(I);
        
    % observed spikes
    ytmp = y0(I);

    % strf prediction
    if ~exist('x','var')
        % fit STRF and make prediction
        nrc = rcSTA(y0(I),X(I,:),length(ops.t),true);
        nrc = normSTA(nrc);
        xtmp = X(I,:) * nrc(:);
        strf(:,:,i) = nrc;
        
    else
        % get prediction and observed spikes
        xtmp = x(I);
    end
    
    % fit each contrast separately
    for j = 1:2
        
        % nonlinearity
        [x_gc{j},y_gc{j},mdl,ahat_gc(j,:),~,n] = ...
            estNL(xtmp(cc==ch(j)),ytmp(cc==ch(j)),ops);
        
        % weighted linear fit
        XX = []; y = [];
        for ii = 1:length(n)
            XX = [XX; [ones(n(ii),1) repmat(x_gc{j}(ii),n(ii),1)]];
            y = [y; repmat(y_gc{j}(ii),n(ii),1)];
        end
        beta_gc(j,:) = regress(y,XX);
        
        % compute slopes and stuff
        xf = linspace(min(x_gc{j}),max(x_gc{j}),100);
        yf = mdl(ahat_gc(j,:),xf);
        maxslope(j) = max(diff(yf)./diff(xf));
        meanslope(j) = mean(diff(yf)./diff(xf));
        
        % prediction in training set
        gc_pred(cc==ch(j)) = mdl(ahat_gc(j,:),xtmp(cc==ch(j)));
        
    end
    
    % threshold predictions by observed max firing rate
    mxfr = max(ytmp)*ops.fs;
    gc_pred(gc_pred>mxfr) = mxfr;
    
    % fit a static model
    [x_static,y_static,mdl,ahat_static,~,n] = estNL(xtmp(:),ytmp(:),ops);
    static_pred = mdl(ahat_static,xtmp);
    static_pred(static_pred>mxfr) = mxfr;
    
    % weighted linear fit (n = number of observations)
    XX = []; y = [];
    for ii = 1:length(n)
        XX = [XX; [ones(n(ii),1) repmat(x_static(ii),n(ii),1)]];
        y = [y; repmat(y_static(ii),n(ii),1)];
    end
    beta_static = regress(y,XX);
    
    % results
    static_ln(i).x = x_static;
    static_ln(i).y = y_static;
    static_ln(i).ahat = ahat_static;
    static_ln(i).pred = static_pred/ops.fs;
    static_ln(i).beta = beta_static;
    gc_ln(i).x = x_gc;
    gc_ln(i).y = y_gc;
    gc_ln(i).ahat = ahat_gc;
    gc_ln(i).maxslope = maxslope;
    gc_ln(i).meanslope = meanslope;
    gc_ln(i).pred = gc_pred/ops.fs;
    gc_ln(i).beta = beta_gc;
    
    fprintf('%d ',i);
        
end
toc;

if exist('x','var');
    % if prediction was supplied, use it
    x_full = x;
else
    % if no prediction was supplied, make one from the mean STRF
    mstrf = mean(strf,3);
    x_full = X * mstrf(:);
end

% gc prediction with average parameters
gc_ahat = mean(cat(3,gc_ln.ahat),3);
gc_pred_full = nan(size(x_full));
xx = cat(1,gc_ln.x); yy = cat(1,gc_ln.y);
for i = 1:2
    gc_pred_full(cc0 == ch(i)) = mdl(gc_ahat(i,:),x_full(cc0==ch(i)))/ops.fs;
    xf = linspace(min([xx{:,i}]),max([xx{:,i}]),100);
    yf = mdl(gc_ahat(i,:),xf);
    gc_mxslope(i) = max(diff(yf)./diff(xf));
    gc_mnslope(i) = mean(diff(yf)./diff(xf));
end    

% static prediction with average parameters
static_ahat = mean(cat(3,static_ln.ahat),3);
static_pred_full = mdl(static_ahat,x_full)/ops.fs;
xf = linspace(min([static_ln.x]),max([static_ln.x]),100);
yf = mdl(static_ahat,xf);
static_mxslope = max(diff(yf)./diff(xf));
static_mnslope = mean(diff(yf)./diff(xf));

% save model object
ops.ln_model = mdl;

% save crossvalidation structs
res.gc_ln_cv = gc_ln;
res.static_ln_cv = static_ln;

% gc model summary
res.gc_ln.x{1} = [xx{:,1}];
res.gc_ln.x{2} = [xx{:,2}];
res.gc_ln.y{1} = [yy{:,1}];
res.gc_ln.y{2} = [yy{:,2}];
res.gc_ln.ahat = gc_ahat;
res.gc_ln.beta = mean(cat(3,gc_ln.beta),3);
res.gc_ln.maxslope = gc_mxslope;
res.gc_ln.meanslope = gc_mnslope;
res.gc_ln.pred = gc_pred_full;
res.gc_ln.corr = corr(gc_pred_full,y0);
res.gc_ln.r2 = res.gc_ln.corr.^2;

% static model summary
res.static_ln.x = [static_ln.x];
res.static_ln.y = [static_ln.y];
res.static_ln.ahat = static_ahat;
res.static_ln.beta = mean([static_ln.beta]',1);
res.static_ln.maxslope = static_mxslope;
res.static_ln.meanslope = static_mnslope;
res.static_ln.pred = static_pred_full;
res.static_ln.corr = corr(static_pred_full,y0);
res.static_ln.r2 = res.static_ln.corr.^2;

% add strf if computed
if exist('strf','var')
    res.static_ln.strf = strf;
end



%ys = mean(reshape(y0,6*ops.fs,[])');
%ps = mean(reshape(gc_pred_full,6*ops.fs,[])');
%figure; hold on;
%plot(ys); plot(ps);