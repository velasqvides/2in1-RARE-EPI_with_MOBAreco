function image = ReconstructImageCartesian(kSpace)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here
bitmask1 = str2double(evalc("bart('bitmask 0 1')"));
image = bart(sprintf('fft -i -u  %i',bitmask1),kSpace);
bitmask1 = str2double(evalc("bart('bitmask 3')"));
image = bart(sprintf('rss %i',bitmask1), image);
end