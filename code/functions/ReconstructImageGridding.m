function image = ReconstructImageGridding(kSpace,traj)
dcf = densityCompRamLak2D(traj);
nSamples = size(kSpace,2);
image = bart(sprintf('nufft -a -d%i:%i:1',nSamples,nSamples), traj, kSpace.*dcf);
bitmask1 = str2double(evalc("bart('bitmask 3')"));
image = bart(sprintf('rss %i',bitmask1), image);
end