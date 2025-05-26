function [imagesEPI_plusRARE,TEs_EPI_plusRARE] = addLastRAREtoEPI_image(imagesRARE,imagesEPI,TEs_EPI,protPara)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ETL_EPI = protPara.ETL_EPI;
TEs_EPI_plusRARE(:,:,:,:,:,2:ETL_EPI + 1) = TEs_EPI;
TEs_EPI_plusRARE(:,:,:,:,:,1) = 0;

nSamples = baseRes * protPara.oversamplingFactor;
imagesEPI_plusRARE = zeros(nSamples,nSamples,1,1,1,ETL_EPI+1,nSlices);
imagesEPI_plusRARE(:,:,:,:,:,2:ETL_EPI + 1,:) = imagesEPI;
imagesEPI_plusRARE(:,:,:,:,:,1,:) = imagesRARE(:,:,:,:,:,end,:);
end