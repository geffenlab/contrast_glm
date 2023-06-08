function r2 = rsquared(obs,pred);

r2 = 1 - sum((obs-pred).^2) ./ sum((obs-mean(obs)).^2);