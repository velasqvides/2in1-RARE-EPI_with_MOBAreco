function [T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE] = reconstruct_MSE_for_hc_004(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE)

[kSpaceMSE,  TEsMSE,  protParaMSE, configMSE] = ...
    preProcessRawDataCartesian(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils);

imagesMSE = ReconstructImageCartesian(kSpaceMSE);
if size(imagesMSE, 6) == 1
    imagesMSE = permute(imagesMSE,[2 1 3 4 5 7 6]);
else
    imagesMSE = permute(imagesMSE,[2 1 3 4 5 6 7]);
end
order = generateReorderingSequence(size(imagesMSE,7));
imagesMSE = imagesMSE(:,:,:,:,:,:,order);
imagesMSE = imagesMSE(:,:,:,:,:,:,2:end-1);
% as(imagesMSE)
iter = 5;
images_SE_abs = abs(imagesMSE); % in case images are reconstructed with phase perservation methods like pics
images_SE_abs_norm = normalizeArray(images_SE_abs);
[M0MSE, R2MSE, T2MSE] = performMobafitRARE(images_SE_abs_norm(:,:,:,:,:,2:end,:), TEsMSE(:,:,:,:,:,2:end), iter);

% binaryMaskRARE=binaryMaskRARE(:,:,2:end-1);
[M0MSE, R2MSE, T2MSE, imagesMSE] = prepareMSEdata(M0MSE, R2MSE, T2MSE,imagesMSE,  binaryMaskRARE, protParaMSE);


end