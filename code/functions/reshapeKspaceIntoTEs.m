function [kSpaceTEs] = reshapeKspaceIntoTEs(kSpace,protPara)
nRARE = protPara.ETL_RARE;
nEPI = protPara.ETL_EPI;
ETL = protPara.ETLfull;
nRAREall = protPara.nSpokesRAREall;
nSpokespPerShot = protPara.nSpokes;
nChannels = protPara.nChannels;
nSlices = protPara.nSlices;
nFE = protPara.baseRes*protPara.oversamplingFactor;
nSpokesFull = nSpokespPerShot*ETL;
% calculate echo line numbers of view ordering
echoViewOrdering = calcEchoViewOrdering(protPara);
% create kSpace data array [nFE,nRViews,nEcho,nCh,nSlices]

kSpaceTEs = zeros(nFE,nSpokespPerShot,nChannels,nSlices,ETL);
for echoNum = 1:ETL
    % use only RARE/EPI echoes with specific TE from each shot
    % retrospective undersampling
    if (protPara.bSameSpokeInET==1)
        if (echoNum <= nRARE)
            % create indices for TE and spoke selection
            subsampIdx = echoViewOrdering(echoNum):(nRARE):nRAREall;
        else
            % create indices for TE and spoke selection
            subsampIdx = nRAREall + echoViewOrdering(echoNum):(nEPI):nSpokesFull;
        end
    else
        % create indices for TE and spoke selection
        subsampIdx = echoViewOrdering(echoNum):(ETL*underSamp):nSpokesFull;
    end

    kSpaceTEs(:,:,:,:,echoNum) = kSpace(:,subsampIdx,:,:);
end
kSpaceTEs = permute(kSpaceTEs,[1 2 3 5 4]);
end
