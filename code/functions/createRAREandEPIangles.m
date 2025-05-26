function [anglesRARE, anglesEPI] = createRAREandEPIangles(protPara)
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
trajAnglesTEs = zeros(nSpokespPerShot,ETL);
for echoNum = 1:ETL
    % use only RARE/EPI echoes with specific TE from each shot
    % retrospective undersampling
    if (protPara.bSameSpokeInET==1)
        if (echoNum <= nRARE)
            % create trajectory
            trajAngles = generateTrajAnglesGDC(nRAREall);
            % create indices for TE and spoke selection
            subsampIdxTraj = echoViewOrdering(echoNum):(nRARE):nRAREall;
        else
            % create trajectory
            trajAngles = generateTrajAnglesGDC(nEPIall);
            % create indices for TE and spoke selection
            subsampIdxTraj = echoViewOrdering(echoNum):(nEPI):nEPIall;
        end
    else
        
        % create indices for TE and spoke selection
        subsampIdxTraj = 1:underSamp:(nSpokesFull/ETL);
    end
    trajAnglesTEs(:,echoNum) = trajAngles(subsampIdxTraj,:);
    
end

ETL_RARE = protPara.ETL_RARE;
% 2. Reorder the trajectroy in a 'BART'-like form
anglesRARE = trajAnglesTEs(:,1:ETL_RARE);
% 3. Reorder the trajectroy in a 'BART'-like form
anglesEPI = trajAnglesTEs(:,(ETL_RARE+1):end);

end
