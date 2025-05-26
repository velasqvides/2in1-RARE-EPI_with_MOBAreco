function [T2maps,M0] = createT2mapsLudger(images,TEs)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nSamples = size(images,1);
nSlices = size(images,7);
T2maps = zeros(nSamples,nSamples,nSlices);
M0 = zeros(nSamples,nSamples,nSlices);
for j =1:nSlices
    [T2maps(:,:,j), M0(:,:,j)] = T2analysis_withM0(abs(images(:,:,1,1,1,:,j)), TEs(:,:,:,:,:,:)*1000);
end

end