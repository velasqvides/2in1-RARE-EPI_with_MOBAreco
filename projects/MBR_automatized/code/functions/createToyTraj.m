function traj = createToyTraj(angles)
ETL = size(angles,2); 
nSamples = 10; 
nSpokes = size(angles,1);
traj = zeros(3,nSamples,nSpokes,1,1,ETL);

for echoN = 1:ETL
    trajAngles = angles(:,echoN);
    writecfl('customAngles',trajAngles);
    traj_ = bart(sprintf('traj -x%i -y%i -r -C customAngles',nSamples,nSpokes));
    traj(:,:,:,:,:,echoN)  = traj_ ; 
end % end for
end % end function