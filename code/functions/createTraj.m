function [traj]= createTraj(angles,ETL,protPara)

nSamples = protPara.baseRes*protPara.oversamplingFactor;
nSpokes = protPara.nSpokes;
traj = zeros(3,nSamples,nSpokes,1,1,ETL);

for echoN = 1:ETL
    trajAngles = angles(:,echoN);
    writecfl('customAngles',trajAngles);
    trajEcho = bart(sprintf('traj -x%i -y%i -r -C customAngles',nSamples,nSpokes));
    traj(:,:,:,:,:,echoN)  = trajEcho; 
end % end for
end % end function