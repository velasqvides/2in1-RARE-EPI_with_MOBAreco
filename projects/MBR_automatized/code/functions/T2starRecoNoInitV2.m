function [T2starOut, B0out] = T2starRecoNoInitV2(kSpace, traj, TEs, protPara, binaryMaskOut, config,specialStamp)

oversamplingFactor = protPara.oversamplingFactor;
baseRes = protPara.baseRes;
nSamples = oversamplingFactor * baseRes;
nChannels = protPara.nChannels;
ETL_EPI = protPara.ETL_EPI;
nSlices = protPara.nSlices;
mobaOversampling = 1.0;
recoSize = baseRes*oversamplingFactor*mobaOversampling; % not sure yet why the factor of 2

bartPath = '~/Documents/tools/bart_v07/bart';
run(fullfile(bartPath, 'startup.m'));


%% 4. compute init (3-point water/fat separation)
mobaOversampling = 1.0;
oversamplingFactor = protPara.oversamplingFactor;
recoSize = baseRes * oversamplingFactor * mobaOversampling; % not sure yet why the factor of 2
R_m0_1feOut = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
R_M1_init_Fout =zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
% recoEPIout = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
% sensEPIout = zeros(recoSize,recoSize,1,nVirtualCoils,nSlices);
waterOut = zeros(baseRes,baseRes,nSlices);
B0out = zeros(baseRes,baseRes,nSlices);
R2starOut = zeros(baseRes,baseRes,nSlices);
T2starOut = zeros(baseRes,baseRes,nSlices);
for slice = 1:nSlices
    TE_e = bart('extract 5 0 3',TEs);
    traj_1fe = bart('extract 5 0 3',traj);
    kdat_1fe = bart('extract 5 0 3',kSpace);
    R_m0_1fe = ...
        bart...
        (sprintf...
        ('moba  -O -G -m0 -i6 -o%4.2f --fat_spec_0 -R2 -t',mobaOversampling),...
        traj_1fe, kdat_1fe(:,:,:,:,:,:,slice), TE_e);
    R_m0_1feOut(:,:,1,1,1,1,:,slice) = R_m0_1fe;

    initB0 = bart('extract 6 2 3',R_m0_1fe);
%     filePath = fullfile(dirToSave,'initR2starOut');
%     initR2starOut = readcfl(filePath);
%     initR2star = initR2starOut(:,:,slice);
    % scale by 0.01 so the inverse will be miliseconds, scale by 1.6 for
    % faster decay
    initR2star = initR2star.*0.01.*1.6;
    initW = lastRAREimageOut(:,:,slice).*1000;
    toPad = (recoSize - baseRes)/2;
    initW = padarray(initW,[toPad toPad]);
    R_M1_init_F = bart('join 6',initW, initR2star, initB0);
    R_M1_init_Fout(:,:,:,:,:,:,:,slice) = R_M1_init_F; 
end

%% 5. moba reconstruction: multi-echo R2* mapping
nSlices=1;
for slice = 1:nSlices
    % moba reconstruction: multi-echo R2* mapping
    [recoEPI, sensEPI] = ...
        bart(...
        sprintf...
        ('moba -G -m3 -rQ:0.1 -rS:0 -rW:3:64:0.1 -rT:3:64:0.1 -u0.01 -o%4.2f  -i10 -C100 -R3 -k --kfilter-2 	 -d4 -t',mobaOversampling),...
        traj, kSpace(:,:,:,:,:,:,slice), TEs);
    recoEPIout(:,:,:,:,:,:,:,slice) = recoEPI;
    sensEPIout(:,:,:,:,slice) = sensEPI;
    cmd1 = sprintf('rm init_file*cfl init_file*hdr');
    system(cmd1);
    %
    tmpMpas = bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes),recoEPI);
    water = bart('slice 6 0',tmpMpas);
    R2star = bart('slice 6 1',tmpMpas);
    B0 = bart('slice 6 2',tmpMpas);

    binaryMask = 1;

    water = binaryMask .* water;
    waterOut(:,:,slice) = water;
    B0 = binaryMask .* B0;
    B0out(:,:,slice) = B0;
    T2star = 1./R2star;
    R2star = binaryMask .* R2star;
    R2starOut(:,:,slice) = R2star;
    T2star = binaryMask .* T2star;
    T2star = bart('scale 1000',T2star);
    % T2star(isinf(T2star)) = 0;
    % T2star(T2star < 0) = 0;
    T2starOut(:,:,slice) = T2star;

    % Create synthesized T2-weighted images
    tmp_result = TEs .* R2star./1000;
    tmp_result1 = bart('scale 1',tmp_result);
    tmp_result = bart('scale  -- -1.0',tmp_result1);
    tmp_exp = bart('zexp',tmp_result);
    synthesizedT2starImages = tmp_exp .* water;
    synthesizedT2starImages = ...
        bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes),synthesizedT2starImages);
    synthesizedT2starImagesOut(:,:,:,:,:,:,slice) = synthesizedT2starImages;

end % end for loop

%% 10. Save the reco data

dirToSave = config.dirToSave;
timeStamp = strrep(strrep(datestr(now),' ','_'),':','-');
finalDirToSave = fullfile(dirToSave,'postprocessed_data','T2star',specialStamp,timeStamp);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave,'recoEPIout');
writecfl(filePath,recoEPIout);
filePath = fullfile(finalDirToSave,'sensEPIout');
writecfl(filePath,sensEPIout);
filePath = fullfile(finalDirToSave,'T2starOut');
writecfl(filePath,T2starOut);
filePath = fullfile(finalDirToSave,'R2starOut');
writecfl(filePath,R2starOut);
filePath = fullfile(finalDirToSave,'waterOut');
writecfl(filePath,waterOut);
filePath = fullfile(finalDirToSave,'B0out');
writecfl(filePath,B0out);
filePath = fullfile(finalDirToSave,'synthesizedT2starImagesOut');
writecfl(filePath,synthesizedT2starImagesOut);
figure, imagesc(flipud(rot90(T2star)), [0 125]);axis equal; axis off; title T2map
end %end function