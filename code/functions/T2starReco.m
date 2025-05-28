function [M0starMaps, R2starMaps, T2starMaps, coilSensEPI, synthesizedEPIimages, B0maps] = T2starReco(kSpace, traj, TEs, initMaps, protPara)

oversamplingFactor = protPara.oversamplingFactor;
baseRes = protPara.baseRes;
nSlices = size(kSpace,7);

mobaOS = protPara.T2StarMobaPara.mobaOversampling;
recoSize = baseRes*oversamplingFactor*mobaOS; 
recoT2starMoba = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
coilSensEPI = zeros(recoSize,recoSize,1,size(kSpace,4),nSlices);
M0starMaps = zeros(recoSize,recoSize,nSlices);
B0maps = zeros(recoSize,recoSize,nSlices);
R2starMaps = zeros(recoSize,recoSize,nSlices);
T2starMaps = zeros(recoSize,recoSize,nSlices);
synthesizedEPIimages = zeros(recoSize,recoSize,1,1,1,size(kSpace,6),nSlices);
binaryMasksEPI = zeros(recoSize,recoSize,nSlices);

regW = protPara.T2StarMobaPara.regW;
regT = protPara.T2StarMobaPara.regT;
u = protPara.T2StarMobaPara.u;
nIterMoba = protPara.T2StarMobaPara.nIterMoba;
nInnerIterMoba = protPara.T2StarMobaPara.nInnerIterMoba;
ReductionFactor = protPara.T2StarMobaPara.ReductionFactor;
B0smoothLevel = protPara.T2StarMobaPara.B0smoothLevel;
B0scaling = protPara.T2StarMobaPara.B0scaling;
sensSmoothLevel = protPara.T2StarMobaPara.sensSmoothLevel;
sensScaling = protPara.T2StarMobaPara.sensScaling;

changeBartVersion(7)
for slice = 1:nSlices
    writecfl('init_file',initMaps(:,:,:,:,:,:,1:3,slice));
    [reco, sens] = ...
        bart(...
        sprintf...
        ('moba -G -m3 -g -rQ:1 -rS:0 -rW:3:64:%f -rT:3:64:%f -u%f -o%f -k --kfilter-2 -i%i -C%i -R%f  -d4 -b%f:%f --sobolev_a %f --sobolev_b %f -I init_file -t',...
        regW, regT, u, mobaOS, nIterMoba, nInnerIterMoba, ReductionFactor, B0smoothLevel, B0scaling, sensSmoothLevel, sensScaling),...
        traj(:,:,:,:,:,:,slice), kSpace(:,:,:,:,:,:,slice), TEs);
    reco = bart(sprintf('resize -c 0 %i 1 %i', recoSize, recoSize), reco);
    sens = bart(sprintf('resize -c 0 %i 1 %i', recoSize, recoSize), sens);
    recoT2starMoba(:,:,:,:,:,:,:,slice) = reco;
    coilSensEPI(:,:,:,:,slice) = sens;
    tmp_maps = squeeze(reco);
    M0 = tmp_maps(:,:,1);
    mask = createBinaryMask(sens, M0, recoSize);
    binaryMasksEPI(:,:,slice) = mask;
    cmd1 = sprintf('rm init_file*cfl init_file*hdr');
    system(cmd1);
    R2star = tmp_maps(:,:,2);
    B0 = tmp_maps(:,:,3);
    M0starMaps(:,:,slice) = M0;
    B0maps(:,:,slice) = B0;
    T2star = 1./R2star;
    R2starMaps(:,:,slice) = R2star;
    T2star = bart('scale 1000',T2star);
    T2starMaps(:,:,slice) = T2star;
    tmp_result = TEs .* R2star./1000;
    tmp_result1 = bart('scale 1',tmp_result);
    tmp_result = bart('scale  -- -1.0',tmp_result1);
    tmp_exp = bart('zexp',tmp_result);
    synthesizedImages = tmp_exp .* M0;
    synthesizedEPIimages(:,:,:,:,:,:,slice) = synthesizedImages;
end 
changeBartVersion(9)
end