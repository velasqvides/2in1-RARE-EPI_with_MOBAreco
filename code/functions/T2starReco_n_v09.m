function [M0starMaps, R2starMaps, T2starMaps, binaryMasksEPI, coilSensEPI, synthesizedEPIimages, B0maps] = T2starReco_n_v09(kSpace, traj, TEs, initMaps, protPara2, ksp_sens, config)

oversamplingFactor = protPara2.oversamplingFactor;
baseRes = protPara2.baseRes;
nSlices = size(kSpace,7);

%% 4. compute init (3-point water/fat separation)
mobaOversampling = protPara2.T2StarMoba.mobaOversampling;
recoSize = baseRes*oversamplingFactor*mobaOversampling*2; % not sure yet why the factor of 2
recoSizeRows = size(initMaps, 1)+4;
recoSizeCols = size(initMaps, 2);
% recoSizeRows = 256;
% recoSizeCols = 256;
recoT2starMoba = zeros(recoSizeRows,recoSizeCols,1,1,1,1,3,nSlices);
coilSensEPI = zeros(recoSizeRows,recoSizeCols,1,size(kSpace,4),nSlices);
M0starMaps = zeros(recoSizeRows,recoSizeCols,nSlices);
B0maps = zeros(recoSizeRows,recoSizeCols,nSlices);
R2starMaps = zeros(recoSizeRows,recoSizeCols,nSlices);
T2starMaps = zeros(recoSizeRows,recoSizeCols,nSlices);
synthesizedEPIimages = zeros(recoSizeRows,recoSizeCols,1,1,1,size(kSpace,6),nSlices);
binaryMasksEPI = zeros(recoSize,recoSize,nSlices);

regW = protPara2.T2StarMobaPara.regW;
regT = protPara2.T2StarMobaPara.regT;
u = protPara2.T2StarMobaPara.u;
nIterMoba = protPara2.T2StarMobaPara.nIterMoba;
nInnerIterMoba = protPara2.T2StarMobaPara.nInnerIterMoba;
ReductionFactor = protPara2.T2StarMobaPara.ReductionFactor;

%% 5. moba reconstruction: multi-echo R2* mapping
changeBartVersion(9)
for slice = 1:nSlices
    % moba reconstruction: multi-echo R2* mapping
    writecfl('init_file',initMaps(:,:,:,:,:,:,:,slice));
    writecfl('ksp_sens_file',ksp_sens);
    [reco, sens] = ...
        bart(...
        sprintf...
        ('moba -G -m3 -g -rQ:1 -rS:0 -rW:3:64:%f -rT:3:64:%f -u%f -o%f -k --kfilter-2 -e0.004 -i%i -C%i -R%f --sobolev_a 110 --sobolev_b 10 -d4 -b10:0.5 -I init_file   -t',regW,regT,u,mobaOversampling,nIterMoba,nInnerIterMoba,ReductionFactor),...
        traj(:,:,:,:,:,:,slice), kSpace(:,:,:,:,:,:,slice), TEs);
    recoT2starMoba(:,:,:,:,:,:,:,slice) = reco;
    coilSensEPI(:,:,:,:,slice) = sens;
    tmp_maps = squeeze(reco);
    % 6. Create mask based on sensitivity maps
    M0 = tmp_maps(:,:,1);
    mask = createBinaryMask(sens, M0, recoSize);
    binaryMasksEPI(:,:,slice) = mask;

    % cmd1 = sprintf('rm init_file*cfl init_file*hdr');
    % system(cmd1);
    %
    R2star = tmp_maps(:,:,2);
    B0 = tmp_maps(:,:,3);

    M0starMaps(:,:,slice) = M0;
    B0maps(:,:,slice) = B0;
    T2star = 1./R2star;
    R2starMaps(:,:,slice) = R2star;
    T2star = bart('scale 1000',T2star);
    % T2star(isinf(T2star)) = 0;
    % T2star(T2star < 0) = 0;
    T2starMaps(:,:,slice) = T2star;

    % Create synthesized T2-weighted images
    tmp_result = TEs .* R2star./1000;
    tmp_result1 = bart('scale 1',tmp_result);
    tmp_result = bart('scale  -- -1.0',tmp_result1);
    tmp_exp = bart('zexp',tmp_result);
    synthesizedImages = tmp_exp .* M0;
%     synthesizedImages = ...
%         bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes),synthesizedImages);
    synthesizedEPIimages(:,:,:,:,:,:,slice) = synthesizedImages;

end % end for loop
% 10. Save the reco data
% savePostProcessedDataEPI(T2starMaps, R2starMaps, M0starMaps, B0maps, binaryMasksEPI, coilSensEPI, recoT2starMoba, synthesizedImages, protPara2, config, initMaps)
% come back to bart v09
changeBartVersion(9)
end %end function