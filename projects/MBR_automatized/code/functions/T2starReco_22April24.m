function [T2starMaps, M0starMaps, R2starMaps, B0maps, binaryMasksEPI] = T2starReco_22April24(kSpace, traj, TEs, protPara, config)

oversamplingFactor = protPara.oversamplingFactor;
baseRes = protPara.baseRes;
nSlices = protPara.nSlices; 

% at the moment T2 MOBA works for bartv07 only
functionDir = fileparts(mfilename('fullpath'));
bartPath = fullfile(functionDir,'../../../../tools/bart_v07/bart');
run(fullfile(bartPath, 'startup.m'));

%% 4. compute init (3-point water/fat separation)
mobaOversampling = protPara.T2StarMoba.mobaOversampling;
recoSize = baseRes*oversamplingFactor*mobaOversampling; % not sure yet why the factor of 2
% % R_m0_1feOut = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
% R_M1_init_Fout =zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
% recoT2starMoba = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
coilSensEPI = zeros(recoSize,recoSize,1,size(kSpace,4),nSlices);
M0starMaps = zeros(recoSize,recoSize,nSlices);
B0maps = zeros(recoSize,recoSize,nSlices);
R2starMaps = zeros(recoSize,recoSize,nSlices);
T2starMaps = zeros(recoSize,recoSize,nSlices);
synthesizedT2starImages = zeros(recoSize,recoSize,1,1,1,size(kSpace,6),nSlices);
binaryMasksEPI = zeros(recoSize,recoSize,nSlices);

regW = protPara.T2StarMobaPara.regW;
regT = protPara.T2StarMobaPara.regT;
u = protPara.T2StarMobaPara.u;
nIterMoba = protPara.T2StarMobaPara.nIterMoba;
nInnerIterMoba = protPara.T2StarMobaPara.nInnerIterMoba;
ReductionFactor = protPara.T2StarMobaPara.ReductionFactor;

%% 5. moba reconstruction: multi-echo R2* mapping
for slice = 1:1
    % moba reconstruction: multi-echo R2* mapping
    
    [reco, sens] = ...
        bart(...
        sprintf...
        ('moba -G -g -m3 -rQ:1 -rS:0 -rW:3:64:%f -rT:3:64:%f -u%f -o%f  -k --kfilter-2 -i%i -C%i -R%f -d4   -t',regW,regT,u,mobaOversampling,nIterMoba,nInnerIterMoba,ReductionFactor),...
        traj, kSpace(:,:,:,:,:,:,slice), TEs);
    recoT2starMoba(:,:,:,:,:,:,:,slice) = reco;
    coilSensEPI(:,:,:,:,slice) = sens;
    tmp_maps = squeeze(reco);
    % 6. Create mask based on sensitivity maps
    M0 = tmp_maps(:,:,1);
    mask = createBinaryMask(sens, M0, recoSize);
    binaryMasksEPI(:,:,slice) = mask;

    cmd1 = sprintf('rm init_file*cfl init_file*hdr');
    system(cmd1);
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
    synthesizedT2starImages(:,:,:,:,:,:,slice) = synthesizedImages;

end % end for loop
% 10. Save the reco data
initMaps=1;
savePostProcessedDataEPI(T2starMaps, R2starMaps, M0starMaps, B0maps, binaryMasksEPI, coilSensEPI, recoT2starMoba, synthesizedImages, protPara, config, initMaps)
% come back to bart v09
functionDir = fileparts(mfilename('fullpath'));
bartPath = fullfile(functionDir,'../../../../tools/bart_v09/bart');
run(fullfile(bartPath, 'startup.m'));
end %end function