function true_gain = compute_true_gain(neuron)

cc = neuron.ops.mean_contrast./neuron.ops.contrast;
cc(isinf(cc)) = 1;

true_gain = zeros(size(cc));
for t = 1:length(cc)
    true_gamma = zeros(size(neuron.ops.beta));
    for lag = 1:size(true_gamma,2)
        lagged_t = t-lag+1;
        if lagged_t<=0
            lagged_t = lagged_t + numel(cc);
        end
        true_gamma(:,lag) = neuron.ops.beta(:,lag) * cc(lagged_t);
    end
    true_gain(t) = 1 + neuron.ops.gain_control * reshape(neuron.ops.beta,[],1)' * ...
        reshape((true_gamma(:,:)-neuron.ops.beta),[],1) / norm(neuron.ops.beta,'fro')^2;
end
