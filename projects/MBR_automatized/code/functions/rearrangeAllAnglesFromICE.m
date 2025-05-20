function [anglesRARE2, anglesEPI2] = rearrangeAllAnglesFromICE(protPara)
nSpokes = protPara.nSpokes;
nRARE = protPara.ETL_RARE;
nEPI = protPara.ETL_EPI;
allAngles = protPara.allAngles;
ETLfull = protPara.ETLfull;
nSlices = protPara.nSlices;

anglesRARE2 = zeros(nSpokes, nRARE);
anglesEPI2 = zeros(nSpokes, nEPI + 1);
for shot = 1:nSpokes
    startIndex = (shot - 1) * ETLfull * nSlices + 1;
    endIndex = startIndex + ETLfull - 1;
    currentAngles = allAngles(startIndex:endIndex);
    
    anglesRARE2(shot, :) = currentAngles(1:nRARE);
    
    anglesEPI2(shot, 1) = currentAngles(nRARE);
   
    anglesEPI2(shot, 2:end) = currentAngles(nRARE+1:end);
end

end