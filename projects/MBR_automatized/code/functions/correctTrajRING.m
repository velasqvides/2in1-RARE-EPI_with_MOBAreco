function [trajCorrected, kSpaceShifts] = correctTrajRING(kSpace,angles)
ETL = size(kSpace,6); 
nSamples = size(kSpace,2); 
nSpokes = size(kSpace,3);
nSlices = size(kSpace,7);
trajCorrected = zeros(3,nSamples,nSpokes,1,1,ETL,nSlices);
kSpaceShifts = zeros(ETL,3,nSlices);
changeBartVersion(9)
for sliceN = 1:nSlices
for echoN = 1:ETL
    tmp = kSpace(:,:,:,:,:,echoN,sliceN); 
    trajAngles = angles(:,echoN);
    writecfl('customAngles',trajAngles);
    traj = bart(sprintf('traj -x%i -y%i -r -C customAngles',nSamples,nSpokes));
%     traj = traj(:,1:255,:);
%     tmp =tmp(:,1:255,:,:);
    spokeShifts = evalc("bart('estdelay -R', traj, tmp)");
    if size(spokeShifts,2) < 32
        spokeShifts = split(spokeShifts,":");
        spokeShifts = arrayfun(@convertCharsToStrings, spokeShifts);
        spokeShifts = arrayfun(@str2num, spokeShifts);
        kSpaceShifts(echoN,:,sliceN) = spokeShifts';
        Sx = spokeShifts(1,1);
        Sy = spokeShifts(2,1);
        Sxy = spokeShifts(3,1);
    else
        spokeShifts = split(spokeShifts,"0m");
        spokeShifts = spokeShifts{2};
        spokeShifts = split(spokeShifts,":");
        spokeShifts = arrayfun(@convertCharsToStrings, spokeShifts);
        spokeShifts = arrayfun(@str2num, spokeShifts);
        kSpaceShifts(echoN,:,sliceN) = spokeShifts';
        Sx = spokeShifts(1,1);
        Sy = spokeShifts(2,1);
        Sxy = spokeShifts(3,1);
    end
    
    tCorrected = bart(sprintf('traj -x%i -y%i -r -O -C customAngles  -q%8.6f:%8.6f:%8.6f',nSamples,nSpokes,Sx,Sy,Sxy));
    
    trajCorrected(:,:,:,:,:,echoN,sliceN)  = tCorrected; 
end % second for
end % first for
end % end function