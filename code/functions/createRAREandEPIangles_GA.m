function [anglesRARE, anglesEPI] = createRAREandEPIangles_GA(protPara)
%viweOrdering: 3.optTinyGA7, 4.tinyGA7, 5.goldenAngles, 6.reorderedGA
nSpokes = protPara.nSpokes;
nRARE = protPara.ETL_RARE;
nEPI = protPara.ETL_EPI;

tau = (sqrt(5) + 1) / 2; % golden ratio
N = 1;
angularSamplingInterval = pi / (tau + N - 1);

index_GA = 0:1:nSpokes - 1;
spokeAngles = angularSamplingInterval * index_GA; % array containing necessary angles for one partition
spokeAngles = mod(spokeAngles, 2 * pi); % projection angles in [0, 2*pi)
angles_RARE = zeros(1, nSpokes * nRARE);

if protPara.ViewOrdering == 3 || protPara.ViewOrdering == 4
    N = 7;
    angularSamplingInterval = pi / (tau + N - 1);
    index_GA = 0:1:nSpokes * nRARE -1;
    angles_RARE = angularSamplingInterval * index_GA;
    angles_RARE = mod(angles_RARE, 2 * pi); % projection angles in [0, 2*pi)
elseif protPara.ViewOrdering == 2
    index_GA = 0:1:nSpokes * nRARE -1;
    angles_RARE = angularSamplingInterval * index_GA;
    angles_RARE = mod(angles_RARE, pi); % projection angles in [0, 2*pi)
elseif protPara.ViewOrdering == 6 %make every second spoke antiparallel
    spokeAngles_pi = mod(spokeAngles, pi);
    [sorted_angles,index_P] = sort(spokeAngles_pi);
    for i = 2:2:length(spokeAngles)
        sorted_angles(i) = mod(sorted_angles(i) + pi, 2*pi);
    end
    for i=1:length(spokeAngles)
        j = index_P(i);
        spokeAngles(j) = sorted_angles(i);
    end
end

if protPara.ViewOrdering == 5 || protPara.ViewOrdering == 6
    index_P = 0:1:nRARE - 1;
    echoRotationAngles = ( (pi / nSpokes) * ((sqrt(5) - 1) / 2) ) * index_P;
    echoRotationAngles = mod(echoRotationAngles, pi/nSpokes);
    selectedPartitions = 1:nRARE;
    partitionIndx = zeros(1, nSpokes * nRARE);
    counter = 1;
    for iR=1:nSpokes
        for iZ=selectedPartitions
            angles_RARE(counter) = spokeAngles(iR) + echoRotationAngles(iZ);
            partitionIndx(counter) = iZ;
            counter = counter + 1;
        end
    end

end

anglesRARE = zeros(nSpokes, nRARE);
for indx = 1 : nRARE
    anglesRARE(:,indx) = angles_RARE(indx:nRARE:end);
end

angles_EPI = zeros(1,nSpokes*nEPI);
if protPara.ViewOrdering == 3
    angleIncrEPI = pi/(nSpokes*(nEPI+1));
else
    angleIncrEPI = pi/(nSpokes*nEPI);
end

anglesLastRAREecho = angles_RARE(nRARE:nRARE:end);
angleIndx = 1;
for shotN = 1:nSpokes
    if protPara.ViewOrdering == 3
        angles_EPI(angleIndx) = anglesLastRAREecho(shotN) + pi - angleIncrEPI;
    else
        angles_EPI(angleIndx) = anglesLastRAREecho(shotN) + pi;
    end
    for iechoN = 1:nEPI-1
        angles_EPI(angleIndx + iechoN) = angles_EPI(angleIndx + iechoN - 1) + pi - angleIncrEPI;
    end
    angleIndx = angleIndx + nEPI;
end

angles_EPI = mod(angles_EPI,2*pi);
anglesEPI = zeros(nSpokes, nEPI + 1);
anglesEPI(:,1) = anglesRARE(:,nRARE);
for indx = 1 : nEPI  
    anglesEPI(:,indx + 1) = angles_EPI(indx:nEPI:end);
end

end % end function
