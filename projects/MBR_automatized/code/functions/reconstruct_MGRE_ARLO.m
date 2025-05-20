function [T2starMGRE, imagesMGRE] = reconstruct_MGRE_ARLO(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE, saveOutput)
[kSpaceMGRE,  TEsMGRE,  protParaMGRE, configMGRE] = ...
    preProcessRawDataCartesian(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils);

imagesMGRE = ReconstructImageCartesian(kSpaceMGRE);
if size(imagesMGRE, 6) == 1
    imagesMGRE = permute(imagesMGRE,[2 1 3 4 5 7 6]);
else
    imagesMGRE = permute(imagesMGRE,[2 1 3 4 5 6 7]);
end

iter = 5;
images_GRE_abs = abs(imagesMGRE);
images_GRE_abs_norm = normalizeArray(images_GRE_abs);
% [M0MGRE, T2starMGRE, R2starMGRE, ~] = performMobafitEPI(images_GRE_abs_norm(:,:,:,:,:,1:end,:), TEsMGRE(:,:,:,:,:,1:end), iter);

images_GRE_abs_norm = squeeze(images_GRE_abs_norm);
TEsMGRE = (squeeze(TEsMGRE))';
for i=1:size(images_GRE_abs_norm,4)
[T2starMGRE(:,:,i), M0MGRE(:,:,i)] = T2analysis_withM0(images_GRE_abs_norm(:,:,:,i), TEsMGRE.*1000);
% figure, imagesc(T2map,[0 150]),axis equal
end
R2starMGRE = 1000./T2starMGRE;
R2starMGRE(isnan(R2starMGRE)) = 0;
T2starMGRE = rearrangeCartesianMaps(T2starMGRE);

[M0MGRE, R2starMGRE, T2starMGRE, imagesMGRE] = prepareMSEdata(M0MGRE, R2starMGRE, T2starMGRE, imagesMGRE, binaryMaskRARE, protParaMGRE);

if saveOutput
    savePostProcessedDataCartesianMGRE(T2starMGRE, R2starMGRE, M0MGRE, protParaMGRE, configMGRE);
    saveAllImagesInPngFormMGRE(T2starMGRE, M0MGRE, imagesMGRE, configMGRE);
end
end