function [T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE] = reconstruct_MGRE_for_hc_004(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE)

[kSpaceMGRE,  TEsMGRE,  protParaMGRE, configMGRE] = ...
    preProcessRawDataCartesian(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils);

imagesMGRE = ReconstructImageCartesian(kSpaceMGRE);
if size(imagesMGRE, 6) == 1
    imagesMGRE = permute(imagesMGRE,[2 1 3 4 5 7 6]);
else
    imagesMGRE = permute(imagesMGRE,[2 1 3 4 5 6 7]);
end
order = generateReorderingSequence(size(imagesMGRE,7));
imagesMGRE = imagesMGRE(:,:,:,:,:,:,order);
imagesMGRE = imagesMGRE(:,:,:,:,:,:,2:end-1);
iter = 5;
images_GRE_abs = abs(imagesMGRE);
images_GRE_abs_norm = normalizeArray(images_GRE_abs);
[M0MGRE, T2starMGRE, R2starMGRE, ~] = performMobafitEPI(images_GRE_abs_norm(:,:,:,:,:,1:end,:), TEsMGRE(:,:,:,:,:,1:end), iter);

[M0MGRE, R2starMGRE, T2starMGRE, imagesMGRE] = prepareMSEdata(M0MGRE, R2starMGRE, T2starMGRE, imagesMGRE, binaryMaskRARE, protParaMGRE);
M0MGRE = rearrangeCartesianMaps(M0MGRE);
R2starMGRE = rearrangeCartesianMaps(R2starMGRE);
T2starMGRE = rearrangeCartesianMaps(T2starMGRE);
imagesMGRE = rearrangeCartesianMaps(imagesMGRE);

end