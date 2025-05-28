function [initMaps] = createInitMpasForT2star(kSpace, traj, TEs, R2, synthesizedRAREimages, binaryMaskRARE, protPara)

mobaOversampling = protPara.T2StarMobaPara.mobaOversampling;
nSamples = size(kSpace, 2)*mobaOversampling;
nSlices = size(kSpace,7);
initMaps = zeros(nSamples,nSamples,1,1,1,1,4,nSlices);

changeBartVersion(7) % at the moment T2 MOBA works better for bartv07 
whichInit = protPara.T2StarMobaPara.whichInit;
scalingFactorM0 = protPara.T2StarMobaPara.scalingFactorM0;
scalingFactorR2star = protPara.T2StarMobaPara.scalingFactorR2star;
mobaOS = protPara.T2StarMobaPara.mobaOversampling;
nIterMobaForInit = protPara.T2StarMobaPara.nIterMobaForInit;
ReductionFactor = protPara.T2StarMobaPara.ReductionFactor;
B0smoothLevel = protPara.T2StarMobaPara.B0smoothLevel;
B0scaling = protPara.T2StarMobaPara.B0scaling;
sensSmoothLevel = protPara.T2StarMobaPara.sensSmoothLevel;
sensScaling = protPara.T2StarMobaPara.sensScaling;

for slice = 1:nSlices
    TE_1    = TEs(:,:,:,:,:,1:3);
    traj_1  = traj(:,:,:,:,:,1:3,slice);
    kdat_1  = kSpace(:,:,:,:,:,1:3,slice);
    init_m0 = bart(sprintf('moba -g -O -G -m0 -o%f -i%i -R%f 2 --fat_spec_0  -b%f:%f --sobolev_a %f --sobolev_b %f -t',...
        mobaOS, nIterMobaForInit, ReductionFactor, B0smoothLevel, B0scaling, sensSmoothLevel, sensScaling),...
        traj_1, kdat_1, TE_1);
    initM0  = bart('extract 6 0 1',init_m0);
    initB0  = bart('extract 6 2 3', init_m0);
    init_m0 = bart(sprintf('moba -g -G -m0 -o%f -i%i -R%f --fat_spec_0  -b%f:%f --sobolev_a %f --sobolev_b %f -t',...
        mobaOS, nIterMobaForInit, ReductionFactor, B0smoothLevel, B0scaling, sensSmoothLevel, sensScaling),...
        traj_1, kdat_1, TE_1);
    initB0_img = bart('extract 6 2 3', init_m0);
    switch whichInit
        case 0
            initM0 = ones(nSamples,nSamples,1);
            initR2star = zeros(nSamples,nSamples,1);
            initB0 = zeros(nSamples,nSamples,1);
            initB0_img = zeros(nSamples,nSamples,1);
        case 1
            lastSynthetic = synthesizedRAREimages(:,:,:,:,:,end,slice).*binaryMaskRARE(:,:,slice);
            initM0 = normalizeArray(lastSynthetic) .* scalingFactorM0;
            initM0 = abs(initM0);
            initR2star = R2(:,:,slice) .* binaryMaskRARE(:,:,slice);
            initR2star = initR2star .* 0.001 .* scalingFactorR2star; % so initR2star is in KiloHertz x scalingFactor
            initB0 = zeros(nSamples,nSamples,1);
            initB0_img = zeros(nSamples,nSamples,1);
        case 2
            initM0 = 0.1*ones(nSamples,nSamples,1);
            initR2star = zeros(nSamples,nSamples,1);
        case 3
            lastSynthetic = synthesizedRAREimages(:,:,:,:,:,end,slice).*binaryMaskRARE(:,:,slice);
            initM0 = normalizeArray(lastSynthetic) .* scalingFactorM0;
            initR2star = R2(:,:,slice).*binaryMaskRARE(:,:,slice);
            initR2star = initR2star .* 0.001 .* scalingFactorR2star; % so initR2star is in KiloHertz x scalingFactor 
        case 4
            lastRAREus = ReconstructImagePics(kSpace(:,:,:,:,:,1,slice), traj(:,:,:,:,:,1,slice)); % first echo in EPI is last RARE
            lastRAREus = lastRAREus .* binaryMaskRARE(:,:,slice);
            lastSynthetic = synthesizedRAREimages(:,:,:,:,:,end,slice) .* binaryMaskRARE(:,:,slice);
            initM0 = sqrt(lastRAREus.^2 + lastSynthetic.^2);
            initM0 = normalizeArray(initM0) .* scalingFactorM0;
            initR2star = R2(:,:,slice).*binaryMaskRARE(:,:,slice);
            initR2star = initR2star .* 0.001 .* scalingFactorR2star;
        case 5
            initR2star = zeros(nSamples,nSamples,1);
    end
    init = bart('join 6',initM0, initR2star, initB0, initB0_img);
    initMaps(:,:,:,:,:,:,:,slice) = init;
end
changeBartVersion(9)
end