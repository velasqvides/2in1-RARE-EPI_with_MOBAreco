function [M0,T2] = T2mappingARLO(img,TEs)
nSlices = size(img,7);
nSamples = size(img,1);
T2 = zeros(nSamples,nSamples,nSlices);
M0 = zeros(nSamples,nSamples,nSlices);
for j = 1:nSlices
[T2(:,:,j), M0(:,:,j)] = T2analysis_withM0(abs(img(:,:,1,1,1,:,j)), TEs.*1000);
end
end