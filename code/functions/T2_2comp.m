function [T2starOut, B0out] = T2_2comp(kSpace, traj, TEs, protPara, config, specialStamp)

oversamplingFactor = protPara.oversamplingFactor;
baseRes = protPara.baseRes;
nSamples = oversamplingFactor * baseRes;
nChannels = protPara.nChannels;
ETL_EPI = protPara.ETL_EPI;
nSlices = protPara.nSlices; 

%% 4. compute init (3-point water/fat separation)
mobaOversampling = 1.0;
recoSize = baseRes*oversamplingFactor*mobaOversampling; % not sure yet why the factor of 2
% R_m0_1feOut = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
R_M1_init_Fout =zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
recoEPIout = zeros(recoSize,recoSize,1,1,1,1,3,nSlices);
sensEPIout = zeros(recoSize,recoSize,1,size(kSpace,4),nSlices);
waterOut = zeros(baseRes,baseRes,nSlices);
B0out = zeros(baseRes,baseRes,nSlices);
R2starOut = zeros(baseRes,baseRes,nSlices);
T2starOut = zeros(baseRes,baseRes,nSlices);
synthesizedT2starImagesOut = zeros(baseRes,baseRes,1,1,1,size(kSpace,6),nSlices);


%% 5. moba reconstruction: multi-echo R2* mapping
for slice = 1:nSlices
    
    [reco, sens] = ...
        bart(...
        sprintf...
        ('moba -G -m2 -rQ:0.1 -rS:0 -rW:3:64:0.1 -u0.01 -o%4.2f  -i8 -C100 -R3 -k --kfilter-2 --fat_spec_0 -d4 -t',mobaOversampling),...
        traj, kSpace(:,:,:,:,:,:,slice), TEs);
    recoOut = reco;
    sensOut = sens;
    

%% 5. Resize outputs
    tmp_maps = bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes),reco);
    tmp_maps = squeeze(tmp_maps);
    tmp_sens = bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes),sens);

    %% 6. Create mask based on sensitivity maps
    bitmask = str2double(evalc("bart('bitmask 3')"));
    tmp_sens_rss = bart(sprintf('rss %i',bitmask),tmp_sens);
    tmp_M0 = tmp_maps(:,:,1);
    tmp_M0_rss_abs = abs(tmp_sens_rss .* tmp_M0);
    circleMask = ...
        createCircleMask([baseRes, baseRes], [baseRes + 1, baseRes + 1]/2, ...
        baseRes/2 - 3);
    tmp_M0_rss_abs = tmp_M0_rss_abs .* circleMask;
    binaryMask = bart('threshold -B 0.5',tmp_M0_rss_abs);
    se = strel('disk',100);
    binaryMask = imclose(binaryMask,se);
    binaryMaskOut(:,:,slice) = binaryMask;

    %
    
    M0_1 = tmp_maps(:,:,1);
    R2_1 = tmp_maps(:,:,2);
    M0_2 = tmp_maps(:,:,3);
    R2_2 = tmp_maps(:,:,4);
    B0   = tmp_maps(:,:,5);

    

    M0_1 = binaryMask .* M0_1;
    M0_2 = binaryMask .* M0_2;
    M0_1Out(:,:,slice) = M0_1;
    M0_2Out(:,:,slice) = M0_2;
    B0 = binaryMask .* B0;
    B0out(:,:,slice) = B0;

    T2_1 = 1./R2_1;
    R2_1 = binaryMask .* R2_1;
    R2_1Out(:,:,slice) = R2_1;
    T2_1 = binaryMask .* T2_1;
    T2_1 = bart('scale 1000',T2_1);
    T2_1Out(:,:,slice) = T2_1;

    T2_2 = 1./R2_2;
    R2_2 = binaryMask .* R2_2;
    R2_2Out(:,:,slice) = R2_2;
    T2_2 = binaryMask .* T2_2;
    T2_2 = bart('scale 1000',T2_2);
    T2_2Out(:,:,slice) = T2_2;

    % Create synthesized T2-weighted images
%     tmp_result = TEs .* R2_2./1000;
%     tmp_result1 = bart('scale 1',tmp_result);
%     tmp_result = bart('scale  -- -1.0',tmp_result1);
%     tmp_exp = bart('zexp',tmp_result);
%     synthesizedT2starImages = tmp_exp .* water;
%     synthesizedT2starImages = ...
%         bart(sprintf('resize -c 0 %i 1 %i',baseRes,baseRes),synthesizedT2starImages);
%     synthesizedT2starImagesOut(:,:,:,:,:,:,slice) = synthesizedT2starImages;

end % end for loop

%% 10. Save the reco data

dirToSave = config.dirToSave;
timeStamp = strrep(strrep(datestr(now),' ','_'),':','-');
finalDirToSave = fullfile(dirToSave,'postprocessed_data','T2_2comp',specialStamp,timeStamp);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave,'recoOut');
writecfl(filePath,recoOut);
filePath = fullfile(finalDirToSave,'sensOut');
writecfl(filePath,sensOut);
filePath = fullfile(finalDirToSave,'T2_1Out');
writecfl(filePath,T2_1Out);
filePath = fullfile(finalDirToSave,'T2_2Out');
writecfl(filePath,T2_2Out);
filePath = fullfile(finalDirToSave,'R2_1Out');
writecfl(filePath,R2_1Out);
filePath = fullfile(finalDirToSave,'R2_2Out');
writecfl(filePath,R2_2Out);
filePath = fullfile(finalDirToSave,'M0_1Out');
writecfl(filePath,M0_1Out);
filePath = fullfile(finalDirToSave,'M0_2Out');
writecfl(filePath,M0_2Out);
filePath = fullfile(finalDirToSave,'B0out');
writecfl(filePath,B0out);

figure, imagesc(flipud(T2_1), [0 150]);axis equal; axis off; title T2_1
figure, imagesc(flipud(T2_2), [0 150]);axis equal; axis off; title T2_2
figure, imagesc(flipud(real(B0)), [-40 40]);axis equal; axis off; title B0
as(M0_1);
as(M0_2);
end %end function