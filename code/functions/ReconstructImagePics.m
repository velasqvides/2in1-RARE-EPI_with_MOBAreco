function image = ReconstructImagePics(kSpace,traj)
% UNTITLED20 Summary of this function goes here
nSamples = size(kSpace,2);
nSlices = size(kSpace,7); % so less slices can be sent to be recostructed from the main script, avoid protPara.nSlices
nEchoes = size(kSpace,6); % same here
image = zeros(nSamples,nSamples,1,1,1,nEchoes,nSlices);
bitmask = str2double(evalc("bart('bitmask 0 1')"));

for j = 1:nSlices
    for i = 1:nEchoes
        coil_img = bart(sprintf('nufft -g -i -d%i:%i:1',nSamples,nSamples), traj(:,:,:,:,:,i), kSpace(:,:,:,:,:,i,j));
        ksp = bart(sprintf('fft -u %i', bitmask), coil_img);
        sens = bart('ecalib -m1 -S', ksp);
        dcf = densityCompRamLak2D(traj);
        writecfl('dcf_file',sqrt(dcf));

        imageEcho = bart(sprintf('pics -g -S -RW:%i:0:0.0075 -RT:%i:0:0.0075 -i50 -e -p dcf_file -t',bitmask,bitmask),traj(:,:,:,:,:,i), kSpace(:,:,:,:,:,i,j), sens);
        % imageEcho = bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes), imageEcho);
        image(:,:,:,:,:,i,j) = imageEcho;
    end
end