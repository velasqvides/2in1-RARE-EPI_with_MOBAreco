function [T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE] = reconstruct_MSE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE)

[kSpaceMSE,  TEsMSE,  protParaMSE, configMSE] = ...
    preProcessRawDataCartesian(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils);

imagesMSE = ReconstructImageCartesian(kSpaceMSE);
if size(imagesMSE, 6) == 1
    imagesMSE = permute(imagesMSE,[2 1 3 4 5 7 6]);
else
    imagesMSE = permute(imagesMSE,[2 1 3 4 5 6 7]);
end
order = generateReorderingSequence(size(imagesMSE,7));
imagesMSE = imagesMSE(:,:,:,:,:,:,order);
if contains(fileName, 'MID00034')
    imagesMSE = imagesMSE(:,:,:,:,:,:,2:end-1);
end
% as(imagesMSE)
iter = 5;
images_SE_abs = abs(imagesMSE); % in case images are reconstructed with phase perservation methods like pics
images_SE_abs_norm = normalizeArray(images_SE_abs);
[M0MSE, R2MSE, T2MSE] = performMobafitRARE(images_SE_abs_norm(:,:,:,:,:,2:end,:), TEsMSE(:,:,:,:,:,2:end), iter);
T2MSE(T2MSE < 0) = 2000;
% M0MSE = rearrangeCartesianMaps(M0MSE);
% R2MSE = rearrangeCartesianMaps(R2MSE);
% T2MSE = rearrangeCartesianMaps(T2MSE);
% imagesMSE = rearrangeCartesianMaps(imagesMSE);
[M0MSE, R2MSE, T2MSE, imagesMSE] = prepareMSEdata(M0MSE, R2MSE, T2MSE,imagesMSE,  binaryMaskRARE, protParaMSE);


end