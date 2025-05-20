function [trajCorrected, kSpaceShifts]= correctTrajRING_playingAround(kSpace,angles,ETL,protPara,gs_shiftsEPI)

nSamples = protPara.baseRes*protPara.oversamplingFactor;
nSpokes = size(kSpace,3);
trajCorrected = zeros(3,nSamples,nSpokes,1,1,ETL);
kSpaceShifts = zeros(ETL,3);
for echoN = 1:ETL
    tmp = kSpace(:,:,:,:,:,echoN); 
    trajAngles = angles(:,echoN);
    writecfl('customAngles',trajAngles);
    traj = bart(sprintf('traj -x%i -y%i -r -c -C customAngles',nSamples,nSpokes));
    % trajTest35 = trajTest34(:,1:255,:);
    % tmp2 =tmp(:,1:255,:,:);
    spokeShifts = evalc("bart('estdelay -R', traj, tmp)");
    if size(spokeShifts,2)<32
        spokeShifts = split(spokeShifts,":");
        spokeShifts = arrayfun(@convertCharsToStrings, spokeShifts);
        spokeShifts = arrayfun(@str2num, spokeShifts);
        kSpaceShifts(echoN,:)=spokeShifts';
        Sx = spokeShifts(1,1);
        Sy = spokeShifts(2,1);
        Sxy = spokeShifts(3,1);
    else
        spokeShifts = split(spokeShifts,"m");
        spokeShifts = spokeShifts{3};
        spokeShifts = split(spokeShifts,":");
        spokeShifts = arrayfun(@convertCharsToStrings, spokeShifts);
        spokeShifts = arrayfun(@str2num, spokeShifts);
        kSpaceShifts(echoN,:)=spokeShifts';
        Sx = spokeShifts(1,1);
        Sy = spokeShifts(2,1);
        Sxy = spokeShifts(3,1);
    end

iseven = rem(ETL, 2) == 0;
if ~iseven
    % writecfl('customAngles',anglesRARE(:,1));
    tCorrected = bart(sprintf('traj -x%i -y%i -r -O -C customAngles  -q%8.6f:%8.6f:%8.6f',nSamples,nSpokes,gs_shiftsEPI(ETL,1)-0.5,gs_shiftsEPI(ETL,2)-0.5,0));
else
    tCorrected = bart(sprintf('traj -x%i -y%i -r -O -C customAngles  -q%8.6f:%8.6f:%8.6f',nSamples,nSpokes,-gs_shiftsEPI(ETL,1),-gs_shiftsEPI(ETL,2),0));
end
    
    trajCorrected(:,:,:,:,:,echoN)  = tCorrected; 
end % end for
end % end function