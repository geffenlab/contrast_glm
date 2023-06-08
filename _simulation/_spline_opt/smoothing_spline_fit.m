function [beta, loss, edf, lam] = smoothing_spline_fit(X, y, S, block, optimize_penalty, lam, gamma, delta)

% add constant term
X = [ones(size(X,1),1), X];
S = blkdiag(1, S);
for j=1:length(block)
    block{j} = block{j}+1;
end
block = [{1}, block];

if nargin < 5
    optimize_penalty = true;
end
if nargin < 6
    lam = repmat(1, length(block), 1);
end
if nargin < 6
    gamma = 1.5;
end
if nargin < 7
    delta = 0.01;
end

maxiter = 1e5;


[beta, mu, loss, z] = pirls_convergence(X, y, S, block, lam, delta);
if ~optimize_penalty
    edf = length(beta);
end

if optimize_penalty
    dGCV = Inf;
    numiter = 0;
    
    tol = 1e-4;
    
    while numiter==0 || (abs(dGCV_old - dGCV) > tol)%(loss_old-loss > eps && numiter < maxiter)
        
        numiter = numiter +1;
        
        loss_old = loss;
        dGCV_old = dGCV;
        
        Slam = S;
        for j=1:length(block)
            Slam(block{j},block{j}) = lam(j)*Slam(block{j},block{j});
        end
        [beta, mu, loss, z] = pirls_step(X, y, Slam, mu);
        
        sqrtW = diag(sqrt(mu));
        [lam, dGCV, edf] = dGCV_convergence(sqrtW*X, sqrtW*z, S, block, gamma, lam);
        
        fprintf("Outer PIRLS loop - iteration %d - dGCV diff %f\n", numiter, dGCV_old - dGCV);
    end
end


end



function [beta, mu, loss, z] = pirls_convergence(X, y, S, block, lam, delta)

maxiter = 1e5;

mu = y + delta;
Slam = S;
for j=1:length(block)
    Slam(block{j},block{j}) = lam(j)*Slam(block{j},block{j});
end

numiter = 0;
loss = Inf;

while numiter==0 || (loss_old-loss > eps && numiter < maxiter)
    
    numiter = numiter + 1;

    loss_old = loss;
    [beta, mu, loss, z] = pirls_step(X, y, Slam, mu);
end

end

function [beta, mu, loss, z] = pirls_step(X, y, Slam, mu)
% prepare pseudodata and weights
eta = log(mu);
z = eta + (y - mu)./mu;
W = diag(mu);
% take step
beta = (X' * W * X + Slam) \ X' * W * z;
% update likelihood (this is only needed to check convergence)
loss = norm(z-X*beta)^2 + beta'*Slam*beta;
% update mu
eta = X * beta;
mu = exp(eta);
end

function [lam, dGCV, edf] = dGCV_convergence(X, y, S, block, gamma, lam)

[Q, R] = qr(X);

tol = 1e-4;

maxiter = 1e5;
numiter = 0;

dGCV = Inf;

while numiter==0 || (dGCV_old-dGCV > tol && numiter < maxiter)
    
    numiter = numiter + 1;

    dGCV_old = dGCV;
    [lam, dGCV, edf] = dGCV_step(Q, R, y, S, block, gamma, lam);
    
    fprintf("Inner dGCV loop - iteration %d - dGCV diff %f\n", numiter, dGCV_old-dGCV);    
end


end

function [lam, dGCV, edf] = dGCV_step(Q, R, y, S, block, gamma, lam)

[dGCV_old, U1, D, V, A] = dGCV_eval(Q, R, y, S, block, gamma, lam);
edf = trace(A);

nlam = length(lam);

y1 = U1' * Q' * y;

n = size(Q,1);
alpha = norm(y-A*y)^2;
delta = n - gamma * edf;


dGCV_grad = zeros(nlam, 1);
dGCV_hess = zeros(nlam);

% dummy values
M_dummy = D \ V' * S * V / D;
F_dummy = M_dummy * (U1' * U1);
M = zeros([size(M_dummy), nlam]);
F = zeros([size(F_dummy), nlam]);

for j = 1:nlam
    lamj = lam(j);
    Sj = zeros(size(S));
    Sj(block{j},block{j}) = S(block{j},block{j});
    M(:,:,j) = D \ V' * Sj * V / D;
    F(:,:,j) = M(:,:,j) * (U1' * U1);
    
    de_delta_de_rho = lamj * trace(F(:,:,j)) * gamma; % note multiplied by gamma
    de_alpha_de_rho = lamj * y1' * (2*M(:,:,j) - F(:,:,j) - F(:,:,j)') * y1;
    
    dGCV_grad(j) = (n/delta^2)*de_alpha_de_rho - (2 * n * alpha / delta^3)*de_delta_de_rho;
end

for j = 1:nlam
    for k = j:nlam
        
        lamk = lam(k);
        Mj = M(:,:,j);
        Mk = M(:,:,k);
        Fj = F(:,:,j);
        Fk = F(:,:,k);
    
        de2_delta_de_rho2 = -2 * lamj * lamk * trace(Mk * Fj) * gamma; % note multiplied by gamma
        
        de2_alpha_de_rho2 = -lamj * lamk * y1' * (2 * Mk * Mj + 2 * Mj * Mk + ...
            - Mj * Fk  - Mk * Fj - Fk' * Mj - Fj' * Mk - Fk * Mj) * y1 +...
            + (k==j) * lamj * y1' * (2 * Mj - Fj' - Fj) * y1;
        
        
        dGCV_hess(j, k) = -(2*n/delta^3)*de_delta_de_rho*de_alpha_de_rho +...
            (n/delta^2)*de2_alpha_de_rho2 - (2*n/delta^3) * de_alpha_de_rho * de_delta_de_rho +...
            (6*n*alpha/delta^4)*de_delta_de_rho*de_delta_de_rho - (2*n*alpha/delta^3)*de2_delta_de_rho2;
        
        dGCV_hess(k, j) = dGCV_hess(j, k);
    end
end
% newton step
rho = log(lam);

[hessV, hessD] = eig(dGCV_hess);
Ds = diag(hessD);
Ds(Ds<=0) = 1e-1;
hessD = diag(Ds);% perturb Hessian to be always positive definite
stepsize = - hessV * hessD * hessV';

rho_test = rho + stepsize \ dGCV_grad;
dGCV = dGCV_eval(Q, R, y, S, block, gamma, exp(rho_test));
while dGCV > dGCV_old
    warning('Newton step was too long. Halving step size...')
    stepsize = stepsize / 2;
    rho_test = rho + stepsize \ dGCV_grad;
    dGCV = dGCV_eval(Q, R, y, S, block, gamma, exp(rho_test));
end
lam = exp(rho_test);

end

function [GCV, U1, D, V, A] = dGCV_eval(Q, R, y, S, block, gamma, lam)

n = size(R,1);

Slam = S;
for j=1:length(block)
    Slam(block{j},block{j}) = lam(j)*Slam(block{j},block{j});
end
Blam = chol(Slam);

[U,D,V] = svd([R; Blam]);

% remove svs that are too small
svs = diag(D);
deficient = svs < max(svs)*sqrt(eps);
D = diag(svs(~deficient));
U = U(:,~deficient);
V = V(:,~deficient);

U1 = U(1:n,:);

A = (Q * U1) * (U1' * Q');
edf = trace(A);

alpha = norm(y-A*y)^2;
delta = n - gamma * edf;

GCV = n * alpha / delta^2;

end