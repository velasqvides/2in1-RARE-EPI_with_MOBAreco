function [kSpaceRARE, kSpaceEPI] = divideRAREandEPIkSpace(kSpaceTEs, protPara)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
%% 1. Divide the RARE data and reorder it in a 'BART-like' form
readoutSamples = protPara.baseRes*protPara.oversamplingFactor;
nSpokes = protPara.nSpokes;
nChannels = protPara.nChannels;
ETL_RARE = protPara.ETL_RARE;
nSlices = protPara.nSlices;
kSpaceRARE = zeros(1,readoutSamples,nSpokes,nChannels,1,ETL_RARE,nSlices);
kSpaceRARE(1,:,:,:,1,:,:) = kSpaceTEs(:,:,:,1:ETL_RARE,:);

%% 2. Divide the EPI data and reorder it in a 'BART-like' form
ETL_EPI = protPara.ETL_EPI;
kSpaceEPI = zeros(1,readoutSamples,nSpokes,nChannels,1,ETL_EPI+1,nSlices);
kSpaceEPI(1,:,:,:,1,:,:) = kSpaceTEs(:,:,:,ETL_RARE:end,:);
end