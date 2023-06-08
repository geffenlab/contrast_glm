function yt = ln_tau(params,tau,ops)

mdl = @(a,x)(a(1) + a(2)*exp((x-a(4))*a(3)));

a = params(:,1);
C = [zeros(1,80) ones(1,80) zeros(1,80)];
t = 1:length(C);
tau = 5;
a_t = a(1) + (a(2) - a(1)).*sum(Ct./length(Ct).*exp(-t./tau))

a_t = a(1) + (a(2) - a(1)).*sum(Ct./length(Ct)).*exp(-t./tau)

n_t = 20;
for t = 1:length(C)
    
    a_t(t) = a(1) + (a(2) - a(1)) .* ...
             sum((C(1:t) ./ t) .* exp(-(1:t)./tau));
end

keyboard

