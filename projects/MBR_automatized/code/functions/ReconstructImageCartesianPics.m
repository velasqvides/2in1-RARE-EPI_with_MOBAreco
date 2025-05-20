function image = ReconstructImageCartesianPics(kSpace)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here
bitmask2 = str2double(evalc("bart('bitmask 0 1')"));
sens = bart('ecalib -W -m1 ', kSpace);
image = bart(sprintf('pics -S -R W:%i:0:0.005 -i40 -e ',bitmask2), kSpace, sens);
% image = bart(sprintf('fft -i -u  %i',bitmask1),kSpace);
bitmask1 = str2double(evalc("bart('bitmask 3')"));
image = bart(sprintf('rss %i',bitmask1), image);
end

%% 256   256     1    12     1    12     7