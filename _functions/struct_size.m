function struct_size(s,scale)

if nargin < 2
    scale = 1;
end
if strcmp(scale,'mb')
    scale = 1e-6;
elseif strcmp(scale,'gb')
    scale = 1e-9;
end

fields = fieldnames(s);

for i = 1:length(fields)
    
    a = s.(fields{i});
    S = whos('a');
    disp({fields{i} S.size S.bytes*scale})
    
end