function p = computeSpikeTrainProb(r,y,bin)

if length(y) ~= length(r)
    error(['Prediction (y) and observed rate (r) must be the same ' ...
           'length!']);
end

p = 1;
for i = 1:length(r)
    
    p = p * ((bin * y(i)) ^ r(i) * exp(-bin*y(i)));
    
end
    
    
    
    