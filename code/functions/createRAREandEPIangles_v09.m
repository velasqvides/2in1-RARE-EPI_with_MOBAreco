function [anglesRARE, anglesEPI] = createRAREandEPIangles_v09(protPara)
%viweOrdering: 3.optTinyGA7, 4.tinyGA7, 5.goldenAngles, 6.reorderedGA
nSpokes = protPara.nSpokes;
nRARE = protPara.ETL_RARE;
nEPI = protPara.ETL_EPI;
tinyGA_RARE = protPara.tinyGA_RARE;
tinyGA_EPI = protPara.tinyGA_EPI;
viewOrdering = protPara.ViewOrdering;
tau = (1 + sqrt(5)) / 2; % golden ratio

if viewOrdering == 1
		angleIncrRARE = pi / (nSpokes * nRARE);
		angleIncrEPI = pi / (nSpokes * (nEPI + 1));
elseif viewOrdering == 2 
		angleIncrRARE = pi / (tau + tinyGA_RARE - 1);
		angleIncrEPI = pi / (tau + tinyGA_EPI - 1);
elseif viewOrdering == 3
		angleIncrRARE = pi / (tau + tinyGA_RARE - 1);
		angleIncrEPI = pi / (nSpokes * (nEPI + 1));	
end

angleIndex = 0:1:(nSpokes * nRARE - 1);
angles_RARE = zeros(1, nSpokes * nRARE);
if viewOrdering == 1
    angles_RARE = angleIncrRARE * angleIndex;
elseif viewOrdering == 2
    angles_RARE = angleIncrRARE * angleIndex;
    angles_RARE = mod(angles_RARE, 2 * pi); % projection angles in [0, 2*pi)
elseif viewOrdering == 3
    angleIndex = 0:1:(nSpokes - 1);
    spokeAngles = angleIncrRARE * angleIndex; % array containing necessary angles for one partition
    spokeAngles = mod(spokeAngles, 2 * pi); % projection angles in [0, 2*pi)
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

anglesLastRAREecho = angles_RARE(nRARE:nRARE:end);

angles_EPI = zeros(1,nSpokes*nEPI);
angleIndx = 1;
for shotN = 1:nSpokes
    angles_EPI(angleIndx) = anglesLastRAREecho(shotN) + pi - angleIncrEPI;
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
