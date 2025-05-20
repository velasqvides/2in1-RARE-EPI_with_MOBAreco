
mainFolder = '/home/cstuser/Users/Velasquez/PhD_project/projects/MBR_automatized/raw_data';
subFolders = dir(mainFolder);
% Filter out the non-directory entries (e.g., '.', '..', and files)
subFolders = subFolders([subFolders.isdir]);
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Exclude '.' and '..'
%%
for i = 7
    folderWitRawData = fullfile(subFolders(i).name);
    filePaths = obtainFilePaths(folderWitRawData);

    %
    numIter = size(filePaths,1);
    for meas = 1:numIter
        fileName = filePaths(meas);
        if contains(fileName, 'RASERHybrid_dev1')
            isOversamplingRemoved = 0;
            nVirtualCoils = 12;
            saveOutput = 1; % set it to 0 when testing
            saveInPng = 1; % set it to 0 when testing

            [kSpaceRARE, kSpaceEPI, trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = ...
                preProcessRawData(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils);
            % imagesRARE = ReconstructImageGridding(kSpaceRARE, trajRARE); % as(imagesRARE)
            % imagesEPI  = ReconstructImageGridding(kSpaceEPI, trajEPI);   % as(imagesEPI)
            % T2 mapping
            protPara.T2MobaPara.mobaOversampling = 1;
            protPara.T2MobaPara.nIterMoba = 10;
            protPara.T2MobaPara.nInnerIterMoba = 250;
            protPara.T2MobaPara.lowerBound = 0.05;
            protPara.T2MobaPara.sensSmoothLevel = 220;
            protPara.T2MobaPara.sensScaling = 15;
            protPara.T2MobaPara.regMoba = 0.0025;

            [M0, R2, T2, binaryMaskRARE, sensRARE, synthesizedRAREimages] = ...
                T2Reco(kSpaceRARE(:,:,:,:,:,2:end,:), trajRARE(:,:,:,:,:,2:end,:), TEsRARE(:,:,:,:,:,2:end), protPara);
            % T2_ = prepareMaps(T2, binaryMaskRARE, protPara);
            % visualizeMaps(T2_,150);
            % visualizeMaps(abs(100.*(T2MSE(:,:,1) - T2)./(T2MSE(:,:,1))),50);

            % T2star mapping
            protPara.T2StarMobaPara.regW = 1;
            protPara.T2StarMobaPara.regT = 1;
            protPara.T2StarMobaPara.u = 0.05;
            protPara.T2StarMobaPara.mobaOversampling = 1;
            protPara.T2StarMobaPara.nIterMoba = 15;
            protPara.T2StarMobaPara.nIterMobaForInit = 6;
            protPara.T2StarMobaPara.nInnerIterMoba = 100;
            protPara.T2StarMobaPara.ReductionFactor = 2.5;
            protPara.T2StarMobaPara.B0smoothLevel = 5;
            protPara.T2StarMobaPara.B0scaling = 0.5;
            protPara.T2StarMobaPara.sensSmoothLevel = 220;
            protPara.T2StarMobaPara.sensScaling = 15;
            protPara.T2StarMobaPara.whichInit = 3;
            protPara.T2StarMobaPara.scalingFactorM0 = 10;
            protPara.T2StarMobaPara.scalingFactorR2star = 2.2;

            initMaps = ...
                createInitMpasForT2star_m3(kSpaceEPI(:,:,:,:,:,1:end,:), trajEPI(:,:,:,:,:,1:end,:), TEsEPI, R2, synthesizedRAREimages, binaryMaskRARE, protPara);

            [M0star, R2star, T2star, sensEPI, synthesizedEPIimages, B0] = ...
                T2starReco_n(kSpaceEPI(:,:,:,:,:,1:end,:), trajEPI(:,:,:,:,:,1:end,:), TEsEPI(:,:,:,:,:,1:end), initMaps, protPara);
            % T2star_ = prepareMaps(T2star, binaryMaskRARE, protPara);
            % visualizeMaps(T2star_,125);
            % visualizeMaps(abs(100.*(T2starMGRE(:,:,1) - T2star)./(T2starMGRE(:,:,1))),50);

            [T2_, R2_, M0_, binaryMaskRARE_, sensRARE_, synthesizedRAREimages_] = ...
                prepareRAREdata(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages,  protPara);
            [T2star_, R2star_, M0star_, sensEPI_, synthesizedEPIimages_, B0_, initB0_, initMaps_] = ...
                prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE, sensEPI, synthesizedEPIimages, B0, initMaps, protPara);

            % visualizeMaps(T2_,150);
            % visualizeMaps(T2star_(:,:,:),125);
            % visualizeMaps(R2star_-10*R2_,20);
            if saveOutput
                savePostProcessedDataRARE(T2_, R2_, M0_, binaryMaskRARE_, sensRARE_, synthesizedRAREimages_, protPara, config);
                savePostProcessedDataEPI(T2star_, R2star_, M0star_, B0_, initB0_, sensEPI_, synthesizedEPIimages_, initMaps_, protPara, config);
            end

            if saveInPng
                saveRAREinPng(T2_, M0_, synthesizedRAREimages_, config);
                saveEPIinPng(T2star_,  M0star_, B0_, initB0_,  synthesizedEPIimages_, binaryMaskRARE_, config);
            end

            %%
        elseif contains(fileName, 'se_mc')
            % meas = 2;
            fileName = filePaths(meas);
            isOversamplingRemoved = 1;
            [T2MSE, imagesMSE] = reconstruct_MSE(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE, saveOutput);
            % visualizeMaps(T2MSE(:,:,:),150);
        elseif contains(fileName, 'gre')
            % meas = 3;
            fileName = filePaths(meas);
            isOversamplingRemoved = 1;
            [T2starMGRE, imagesMGRE] = reconstruct_MGRE(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE, saveOutput);
            % visualizeMaps(T2starMGRE(:,:,:),125);
        else
            % If file doesn't match any of the conditions, you can skip or handle it differently
            disp(['Skipping unknown file: ', fileName]);
        end
    end
end
