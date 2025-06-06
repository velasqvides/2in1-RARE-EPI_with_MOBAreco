function [circleMask] = createCircleMask(dim, center, radius)
% function created by Ludger Starke
[x,y] = meshgrid(-(center(2)-1):(dim(2)-center(2)),-(center(1)-1):(dim(1)-center(1)));
circleMask = ((x.^2+y.^2) <= radius^2);

end