 function [trajRARE, trajEPI] = createRAREandEPItrajectories(protPara)
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here
nRARE = protPara.ETL_RARE;
nEPI = protPara.ETL_EPI;
ETL = protPara.ETLfull;
nRAREall = protPara.nSpokesRAREall;
nEPIall = protPara.nSpokesEPIall;
nSpokespPerShot = protPara.nSpokes;
nFE = protPara.baseRes * protPara.oversamplingFactor;
nSpokesFull = nSpokespPerShot*ETL;
% calculate echo line numbers of view ordering
echoViewOrdering = calcEchoViewOrdering(protPara);
% create kSpace data array [nFE,nRViews,nEcho,nCh,nSlices]
trajTEs = zeros(3,nFE,nSpokespPerShot,1,1,ETL);
for echoNum = 1:ETL
    % use only RARE/EPI echoes with specific TE from each shot
    % retrospective undersampling
    if (protPara.bSameSpokeInET==1)
        if (echoNum <= nRARE)
            % create trajectory
            kSpaceTraj = generatekSpaceTrajectory(nRAREall,nFE);
            % create indices for TE and spoke selection
            subsampIdxTraj = echoViewOrdering(echoNum):(nRARE):nRAREall;
        else
            % create trajectory
            kSpaceTraj = generatekSpaceTrajectory(nEPIall,nFE);
            % create indices for TE and spoke selection
            subsampIdxTraj = echoViewOrdering(echoNum):(nEPI):nEPIall;
        end
    else
        % create trajectory
        kSpaceTraj = generatekSpaceTrajectory(nSpokesFull/ETL,nFE);
        % create indices for TE and spoke selection
        subsampIdxTraj = 1:underSamp:(nSpokesFull/ETL);
    end

    trajTEs(:,:,:,1,1,echoNum) = kSpaceTraj(:,:,subsampIdxTraj);
end

ETL_RARE = protPara.ETL_RARE;
% 2. Reorder the trajectroy in a 'BART'-like form
trajRARE = trajTEs(:,:,:,1,1,1:ETL_RARE);
% 3. Reorder the trajectroy in a 'BART'-like form
trajEPI = trajTEs(:,:,:,1,1,(ETL_RARE+1):end);

end
