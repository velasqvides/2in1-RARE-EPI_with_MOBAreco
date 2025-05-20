function [kSpaceTEs] = reshapeKspaceIntoTEs_GA(kSpace,protPara)
ETL = protPara.ETLfull;
nSpokespPerShot = protPara.nSpokes;
nChannels = protPara.nChannels;
nSlices = protPara.nSlices;
nFE = protPara.baseRes*protPara.oversamplingFactor;
nSpokesFull = nSpokespPerShot*ETL;
% calculate echo line numbers of view ordering
% create kSpace data array [nFE,nRViews,nEcho,nCh,nSlices]
% echoViewOrdering=1:32;
kSpaceTEs = zeros(nFE,nSpokespPerShot,nChannels,nSlices,ETL);
for echoNum = 1:ETL
    subsampIdx = echoNum:(ETL):nSpokesFull;
    kSpaceTEs(:,:,:,:,echoNum) = kSpace(:,subsampIdx,:,:);
end
kSpaceTEs = permute(kSpaceTEs,[1 2 3 5 4]);
end
