function [T2MSE, imagesMSE] = reconstruct_MSE_ARLO(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE, saveOutput)

[kSpaceMSE,  TEsMSE,  protParaMSE, configMSE] = ...
    preProcessRawDataCartesian(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils);


imagesMSE = ReconstructImageCartesian(kSpaceMSE);
if size(imagesMSE, 6) == 1
    imagesMSE = permute(imagesMSE,[2 1 3 4 5 7 6]);
else
    imagesMSE = permute(imagesMSE,[2 1 3 4 5 6 7]);
end
% as(imagesMSE)
iter = 5;
images_SE_abs = abs(imagesMSE); % in case images are reconstructed with phase perservation methods like pics
images_SE_abs_norm = normalizeArray(images_SE_abs);
images_SE_abs_norm = squeeze(images_SE_abs_norm);
TEsMSE = (squeeze(TEsMSE))';
for i=1:size(images_SE_abs_norm, 4)
[T2MSE(:,:,i), M0MSE(:,:,i)] = T2analysis_withM0(images_SE_abs_norm(:,:,2:end,i), TEsMSE(1,2:end).*1000);
% figure, imagesc(T2map,[0 150]),axis equal
end
R2MSE = 1000./T2MSE;
R2MSE(isnan(R2MSE)) = 0;
T2MSE = rearrangeCartesianMaps(T2MSE);
% [M0MSE, R2MSE, T2MSE] = performMobafitRARE(images_SE_abs_norm(:,:,:,:,:,2:end,:), TEsMSE(:,:,:,:,:,2:end), iter);
% T2MSE = rearrangeCartesianMaps(T2MSE);

[M0MSE, R2MSE, T2MSE, imagesMSE] = prepareMSEdata(M0MSE, R2MSE, T2MSE,imagesMSE,  binaryMaskRARE, protParaMSE);

if saveOutput
    savePostProcessedDataCartesian(T2MSE, R2MSE, M0MSE, protParaMSE, configMSE);
    saveAllImagesInPngFormMSE(T2MSE, M0MSE, imagesMSE, configMSE);
end
end