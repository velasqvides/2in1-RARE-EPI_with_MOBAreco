%gardient delay exmaples with bar trajectories
delayX = -0.6;
delayY = -0.9;
nSpokes=20;
nSamples = 256;
tnomi = bart(sprintf('traj -x%i -y%i -r',nSamples,nSpokes));
tCorrupted = bart(sprintf('traj -x%i -y%i -r -O -q %f:%f:0',nSamples,nSpokes,delayX,delayY));
angles=0:pi/nSpokes:(pi/nSpokes)*(nSpokes-1);
% angles=0:pi/4:(3*pi/4);

k_shifts = (cos(angles).^2).*(delayX) + (sin(angles).^2).*(delayY) ;
figure,quiver2Dspokes(tnomi)
figure,quiver2Dspokes(tCorrupted)
corrected = zeros(3,nSamples,nSpokes);
for i=1:nSpokes
corrupted = tCorrupted(:,:,i);
corrected(:,:,i)=corrupted+[sin(angles(i));cos(angles(i));0].*[-delayY;-delayX;0];
end
figure,quiver2Dspokes(corrected)