function g = gain(obj,contrast,ops)

g = obj.gain_control * obj.mean_contrast./contrast +...
    (1-obj.gain_control) * 1;
g(isinf(g)) = 1;
g(isnan(g)) = 1;

if isfield(obj,'tau') && ~isempty(obj.tau) && all(obj.tau > 0) && obj.gain_control ~= 0
    
    % get gain values
    gg = unique(g,'stable');
    
    % get second contrast
    ntt = sum(g==gg(1)) / ( ops.exemplars*ops.ntrials );
    x = 1:ntt;
    
    % define decay from second to first contrast
    % y = end + (start - end)e^(-tau*x)
    g_new = obj.expfun([gg(2),gg(1),obj.tau(1)],x);
    g(g == gg(1)) = repmat(g_new,1,ops.exemplars*ops.ntrials);
    
    if numel(obj.tau) == 2
        
        % add transition from first to second contrast
        ntt = sum(g==gg(2)) / ( ops.exemplars*ops.ntrials );
        x = 1:ntt;
        
        % define decay
        % y = end + (start - end)e^(-tau*x)
        g_new = obj.expfun([gg(1),gg(2),obj.tau(2)],x);
        g(g == gg(2)) = repmat(g_new,1,ops.exemplars*ops.ntrials);
        
    end

end

