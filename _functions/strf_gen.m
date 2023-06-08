function beta = strf_gen(ops)
% model strf as 2d gaussian
mu = [20, 2];
sig = [0.8 0.1; 0.1 .5];
[X1,X2] = meshgrid(1:length(ops.f),1:length(ops.t));
X = [X1(:) X2(:)];
yy = mvnpdf(X,mu,sig);
realSTA = reshape(yy,length(ops.t),length(ops.f))';
noise = .025*randn(size(realSTA,1),size(realSTA,2));
%noise = 1*randn(size(realSTA,1),size(realSTA,2));
realSTA = realSTA + noise;
realSTA = realSTA - mean(realSTA(:));
beta = realSTA ./ sqrt(sum(realSTA(:).^2));