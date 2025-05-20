function [M0maps, R2maps, T2maps, binaryMaskRARE, sensRARE, synthesizedRAREimages] = T2Reco(kSpace, traj, TEs, protPara)
nSlices = size(kSpace,7);
oversamplingFactor = protPara.oversamplingFactor;
baseRes = protPara.baseRes;
ETL_RARE = size(kSpace,6);
mobaOversampling = protPara.T2MobaPara.mobaOversampling;
recoSize = baseRes*oversamplingFactor*mobaOversampling; % not sure yet why the factor of 2
recoT2Moba = zeros(recoSize,recoSize,1,1,1,1,2,nSlices);
sensRARE = zeros(recoSize,recoSize,1,size(kSpace,4),nSlices);
T2maps = zeros(recoSize,recoSize,nSlices);
R2maps = zeros(recoSize,recoSize,nSlices);
M0maps = zeros(recoSize,recoSize,nSlices);
synthesizedRAREimages = zeros(recoSize,recoSize,1,1,1,ETL_RARE,nSlices);
binaryMaskRARE = zeros(recoSize,recoSize,nSlices);

nIterMoba = protPara.T2MobaPara.nIterMoba;
nInnerIterMoba = protPara.T2MobaPara.nInnerIterMoba;
regMoba = protPara.T2MobaPara.regMoba;
lowerBound = protPara.T2MobaPara.lowerBound;
sensSmoothLevel = protPara.T2MobaPara.sensSmoothLevel;
sensScaling = protPara.T2MobaPara.sensScaling;

changeBartVersion(7) % at the moment T2 MOBA works for bartv07 only

for slice = 1:nSlices
    [reco, sens] = ...
        bart(sprintf('moba -F -g -i%i -n -C%i -j%f -o%f -d4  -B%f -k --kfilter-2 --sobolev_a %f --sobolev_b %f -t',...
        nIterMoba, nInnerIterMoba, regMoba, mobaOversampling, lowerBound, sensSmoothLevel, sensScaling), ...
        traj(:,:,:,:,:,:,slice), kSpace(:,:,:,:,:,:,slice), TEs);

    recoT2Moba(:,:,:,:,:,:,:,slice) = reco;
    sensRARE(:,:,:,:,slice) = sens;
    tmp_maps = squeeze(reco);
    % 6. Create mask based on sensitivity maps
    M0 = tmp_maps(:,:,1);
    mask = createBinaryMask(sens, M0,recoSize);
    binaryMaskRARE(:,:,slice) = mask;
    % 8. M0, R2 and T2
    M0maps(:,:,slice) = M0;
    R2 = tmp_maps(:,:,2);
    R2 = R2 .* 10; % so R2 is in Herz now
    R2maps(:,:,slice) = R2;
    T2 = 1./R2;
    T2 = bart('scale 1000',T2); % so T2 is in ms
    T2maps(:,:,slice) = T2;
    % 9. Create synthesized T2-weighted images
    tmp_result = TEs .* R2;
    tmp_result = bart('scale  -- -1.0',tmp_result);
    tmp_exp = bart('zexp',tmp_result);
    SynthesizedImages = tmp_exp .* M0;
    synthesizedRAREimages(:,:,:,:,:,:,slice) = SynthesizedImages;
end
changeBartVersion(9)
end