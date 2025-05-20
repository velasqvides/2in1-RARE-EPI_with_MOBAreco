function image = ReconstructImageNufft(kSpace,traj)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here

nSamples = size(kSpace,2);
image = bart(sprintf('nufft -i -d%i:%i:1',nSamples,nSamples), traj, kSpace);
bitmask1 = str2double(evalc("bart('bitmask 3')"));
image = bart(sprintf('rss %i',bitmask1), image);
end