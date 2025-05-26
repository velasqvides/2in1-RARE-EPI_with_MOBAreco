% generate trajectory for reconstruction
function [kSpaceTrajectory] = generatekSpaceTrajectory(nSpokesFull,baseResolution)

% kSpaceTrajectory = bart(sprintf('traj -x%i -y%i -r',baseResolution,nSpokesFull));
kSpaceTrajectory = bart(sprintf('traj -x%i -y%i -r', baseResolution,nSpokesFull));
end
