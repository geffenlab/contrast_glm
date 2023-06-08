clear all; close all;
chrisComp = false;
[~,name] = system('hostname');
if contains(name,'oto-gefen-007.MED.UPENN.EDU')
    chrisComp = true;
    addpath(genpath('~/chris-lab/code_general/'));
end
%addpath(genpath('/Users/chris/chris-lab/code_general/fitting/glmnet_matlab/'));



%%%%%%%%%%
%% options
ops.w = .5;
ops.period = .025;
ops.fs = 1/ops.period;
ops.f = 4000*(2 .^ ([0:32]/8));
ops.t = 0:1/ops.fs:(ops.w-ops.period);
ops.mu = [0 0];
ops.sd = [1 3];
ops.exemplars = 20;
ops.blockLength = 2;
ops.ntrials = 20;
ops.pad = ops.w;



%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the stimulus
[stim,c] = stimgen(ops);
contrast = c; %std(stim,1,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% forward model simulation
obj.base_rate = 10;
obj.contrast = contrast;
obj.mean_contrast = 1/((ops.sd(1)+ops.sd(2))/(2*ops.sd(1)*ops.sd(2)));
obj.gain_control = .5;
obj.operating_point = ops.mu(1);
obj.beta = strf_gen(ops);

% run simulation
[y,l,h,g,x0] = simulate(obj,stim,contrast);
mx = mean(stim,2);



%%%%%%%%%%%%%%
%% fit the glm
% make the stim predictors
x = lagDesignMatrix(stim-obj.operating_point,length(ops.t))';

% make the contrast predictor
cc = obj.mean_contrast./contrast;
cc(isinf(cc)) = 1;

% make the interaction predictor
xc = lagDesignMatrix((stim-obj.operating_point).*cc,length(ops.t))';

% combine
dm = [x,xc,cc];

% parallel pool for fitting
delete(gcp('nocreate'));
p = parpool('local',3);

% fit
fprintf('glmnet fit... '); tic;
opts = glmnetSet;
opts.alpha = .95; % elastic net, close to lasso
fit = cvglmnet(dm,y,'poisson',opts,[],[],[],true); toc;



%%%%%%%%%%%%%%%%%%
%% post processing

% extract fit coefficients and prediction
coeffs = cvglmnetPredict(fit,[],'lambda_min','coefficients');
pred = cvglmnetPredict(fit,dm,'lambda_min','response');
beta0 = coeffs(1);
beta1 = coeffs(2:2+length(ops.f)*length(ops.t)-1);
beta2 = coeffs(2+length(ops.f)*length(ops.t):length(coeffs)-1);
beta3 = coeffs(end);

% reshape for visualization
b1 = reshape(beta1,length(ops.f),[]);
b2 = reshape(beta2,length(ops.f),[]);

% contrast design matrix
C = lagDesignMatrix(repmat(cc,1,length(ops.f)),length(ops.t));

% gain index (not sure if this is correct)
% w = 1 + beta2'./(beta1' + beta2') * (C-1);

% correct averaging (matrix sum, normalized to appropriately
% average):
bs = beta1 + beta2;
bs(bs==0) = 1e-5;
w = 1 + (beta2'./bs'*(C-1)) / ...
     (length(ops.f)*length(ops.t));
wm = reshape(w',(ops.blockLength*2 + ops.pad)*ops.fs,[])';



%%%%%%%%%%%
%% plotting
f1 = figure(1); clf
subplot(3,4,1); hold on;
imagesc(ops.t,ops.f/1000,obj.beta); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('STRF (beta)'); axis tight;

subplot(3,4,2); hold on;
scatter(h(c==ops.sd(2)),y(c==ops.sd(2)),'r.');
scatter(h(c==ops.sd(1)),y(c==ops.sd(1)),'b.');
xlabel('Linear Drive');
ylabel('Firing Rate');
title('Nonlinearity')

subplot(3,4,3); hold on;
imagesc(reshape(h,180,[])'); colorbar;
plot(repmat(ops.blockLength*ops.fs,2,1),[0 ops.ntrials*ops.exemplars],'r','linewidth',2)
xlabel('Time (step)');
ylabel('Trial');
title(sprintf('Linear Drive\nc_{low}=%02.1f,c_{hi}=%02.1f',ops.sd(1),ops.sd(2)));
axis tight;

subplot(3,4,4); hold on;
imagesc(reshape(y,180,[])'); colorbar;
plot(repmat(ops.blockLength*ops.fs,2,1),[0 ops.ntrials*ops.exemplars],'r','linewidth',2)
plot(mean(reshape(y,180,[])'),'w','linewidth',1)
xlabel('Time (step)');
ylabel('Trial');
title('Poisson Counts');
axis tight;

subplot(3,4,5); hold on;
imagesc(ops.t,ops.f/1000,b1); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('${\beta_1}$','Interpreter','Latex'); axis tight;

subplot(3,4,6); hold on;
imagesc(ops.t,ops.f/1000,b2); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('${\beta_2}$','Interpreter','Latex'); axis tight;

subplot(3,4,7); hold on;
imagesc(ops.t,ops.f/1000,b1+b2); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('${\beta}_1 + {\beta}_2$',...
      'Interpreter','latex'); axis tight;

subplot(3,4,8); hold on;
imagesc(ops.t,ops.f/1000,b2./(b1+b2)); colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title('${{\beta}_2 / ({\beta}_1 + {\beta}_2)}$',...
      'Interpreter','latex'); axis tight;

subplot(3,4,9); hold on;
plot(mean(wm),'k');
plot(g(1:180),'k--');
ylabel('$\bar{w},g$','Interpreter','Latex');
title('Gain Modulation')
lh = legend('$\bar{w} = 1 + ({\beta}_2 / ({\beta}_1 + {\beta}_2) \cdot (C-1)) / (T*F)$',...
            '$g$',...
            'Interpreter','latex');
xlabel('Time (step)');
lh.Position(1) = lh.Position(1) + .2;


if chrisComp
    saveFigPDF(f1,[1000 600],...
               sprintf('sim_res_g%de-1_fr%d_c%d-%de-1.pdf',...
                       obj.gain_control*10, ...
                       obj.base_rate, ...
                       ops.sd(1)*10,ops.sd(2)*10));
end























%% functions


% gain function
function g = gain(obj,contrast)
    g = obj.gain_control * obj.mean_contrast./contrast +...
        (1-obj.gain_control) * 1;
    g(isinf(g)) = 1;
end



% gain scaling function
function [lambda_scaled, this_gain] = gain_scaling(obj,lambda,contrast)
    this_gain = gain(obj,contrast);
    lambda_scaled = log(obj.base_rate) + this_gain .* (lambda-log(obj.base_rate));
end



% firing rate function
function [l,linear_drive,g,x0] = lambda(obj,x,contrast)
    % x = design matrix
    % beta = strf weights    
    x0 = convSTRF(x'-obj.operating_point,fliplr(obj.beta))';
    linear_drive = log(obj.base_rate) + x0;
    [drive_with_gain,g] = gain_scaling(obj,linear_drive,contrast);
    l = exp(drive_with_gain);
    % l(isinf(l)) = 0;
end



% simulation
function [y,l,h,g,x0] = simulate(obj,x,contrast)
    % y = spike count
    % l = spike rate
    % h = linear drive
    % g = gain
    [l,h,g,x0] = lambda(obj,x,contrast);
    y = poissrnd(l);
end



% generate a model strf
function beta = strf_gen(ops)
    % model strf as 2d gaussian
    mu = [20, 2];
    sig = [0.8 0.1; 0.1 .5];
    [X1,X2] = meshgrid(1:length(ops.f),1:length(ops.t));
    X = [X1(:) X2(:)];
    yy = mvnpdf(X,mu,sig);
    realSTA = reshape(yy,length(ops.t),length(ops.f))';
    noise = .025*randn(size(realSTA,1),size(realSTA,2));
    realSTA = realSTA + noise;
    realSTA = realSTA - mean(realSTA(:));
    beta = realSTA ./ sqrt(sum(realSTA(:).^2));
end



% compute gain modulation index
function w = gain_modulation_index(obj,beta1,beta2,contrast)
    % beta 1 = stimulus weights (STRF?)
    % beta 2 = contrast interaction weights (gain field?)
    w = 1 + beta2 / (beta1+beta2) .* ...
        (obj.mean_contrast./contrast-1);
end



% 1D convolution for multiple bands
function convSS = convSTRF(S,STA,mode)

    % function convSS = convSTRF(S,STA)
    % computes the linear combination of STA (neural filter) and the
    % stimulus (S)
    if ~exist('mode','var')
        mode = 'full';
    end
    convSS = zeros(1,size(S,2));
    % for each frequency
    for i = 1:size(S,1)
        % compute convolution of corresponding STA band and stim band
        x = conv(S(i,:), fliplr(STA(i,:)), mode);
        %x = x(floor(length(STA)/2)+1:end-floor(length(STA)/2));
        if strcmp(mode,'full')
            x = x(1:size(S,2));
            % account for conv shift by 1 sample
            x = [0 x(1:end-1)];
        end
        convSS = convSS + x;
    end
end



% design matrix lag
function X = lagDesignMatrix(S,lags)

    if size(S,2) < size(S,1)
        S = S';
    end

    nfs = size(S,1);
    X = zeros(nfs*lags,size(S,2));
    for i = 1:lags
        
        rowI = (i-1)*nfs+1:i*nfs;
        X(rowI,:) = circshift(S,(i-1),2);

    end

end



% figure code
function saveFigPDF(h,sz,savepath);
    if size(sz,1) == 2
        sz = sz';
    end

    set(groot,{'DefaultAxesXColor','DefaultAxesYColor', ...
               'DefaultAxesZColor'},{'k','k','k'})
    set(h,'PaperPositionMode','auto');         
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','points');
    if ~isempty(sz)
        set(h,'PaperSize',sz);
        set(h,'Position',[0 0 sz]);
    end
    if exist('savepath','var') & ~isempty(savepath)
        print(h,savepath,'-dpdf','-r300','-painters');
    end
end

    
