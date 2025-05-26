function [kSpaceEPI_plusRARE,trajEPI_plusRARE,TEs_EPI_plusRARE] = addLastRAREtoEPI_data(kSpaceRARE,kSpaceEPI,trajRARE,trajEPI,TEs_EPI,protPara)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ETL_EPI = protPara.ETL_EPI;
TEs_EPI_plusRARE(:,:,:,:,:,2:ETL_EPI + 1) = TEs_EPI;
TEs_EPI_plusRARE(:,:,:,:,:,1) = 0;
nSlices = size(kSpaceEPI,7); % this function needs to serve even when only one slice is sent 
nSamples = protPara.baseRes * protPara.oversamplingFactor;
kSpaceEPI_plusRARE = zeros(1,nSamples,protPara.nSpokes,20,1,ETL_EPI+1,nSlices);
kSpaceEPI_plusRARE(:,:,:,:,1,2:ETL_EPI + 1,:) = kSpaceEPI;
kSpaceEPI_plusRARE(:,:,:,:,:,1,:) = kSpaceRARE(:,:,:,:,:,end,:);

trajEPI_plusRARE = zeros(3,nSamples,protPara.nSpokes,1,1,ETL_EPI+1);
trajEPI_plusRARE(:,:,:,1,1,2:ETL_EPI + 1) = trajEPI;
trajEPI_plusRARE(:,:,:,1,1,1) = trajRARE(:,:,:,:,:,end);
end