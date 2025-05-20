function [ mask ] = createBackgroundMask( dim,side )

mask = false(dim);
mask(1:side,1:side) = true;
mask(end-side:end,1:side) = true;
mask(1:side,end-side:end) = true;
mask(end-side:end,end-side:end) = true;

end

