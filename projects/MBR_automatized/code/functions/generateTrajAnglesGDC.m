% generate trajectory for reconstruction
function [trajAngles] = generateTrajAnglesGDC(nSpokesFull)
% define angle increment between spokes
angleIncr=pi/nSpokesFull;
trajAngles = zeros(nSpokesFull,1);
for i=1:nSpokesFull
    trajAngles(i,1) = angleIncr*(i-1);  
end

end
