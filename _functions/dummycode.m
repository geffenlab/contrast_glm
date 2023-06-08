function [C C1] = dummycode(C,cc,ops)

C1 = C;
if ops.dummycode > 0
    
    % ops.shift = [];
    if ~isempty(ops.shift)
        % shift to cover transition time
        cs = circshift(cc,ops.shift);
    else
        cs = cc;
    end
    
    % split contrast predictors
    uc = unique(cc,'stable');
    c1 = C; c2 = C;
    c1(cs == uc(2),:) = 0;
    c2(cs == uc(1),:) = 0;
    
    if ops.dummycode == 2

        % dummy code beta3
        C1 = [c1 c2];
        
    end
    C = [c1 c2];
    
end