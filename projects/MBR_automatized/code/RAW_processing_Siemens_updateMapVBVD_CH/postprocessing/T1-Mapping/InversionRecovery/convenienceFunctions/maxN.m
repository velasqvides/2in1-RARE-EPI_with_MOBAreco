function out = maxN(in)
% Generalization of the sum function to an arbitrary number of dimensions

dim = numel(size(in));

out = in;

for ii = 1:dim
    
    out = squeeze(max(out));
    
end