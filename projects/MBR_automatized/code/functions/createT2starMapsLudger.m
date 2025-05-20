function [T2starMaps,M0star] = createT2starMapsLudger(images,TEs)

nSamples = size(images,1);
nSlices = size(images,7);
T2starMaps = zeros(nSamples,nSamples,nSlices);
M0star = zeros(nSamples,nSamples,nSlices);

for j =1:nSlices
    [T2starMaps(:,:,j), M0star(:,:,j)] = T2analysis_withM0(abs(images(:,:,1,1,1,1:end,j)), TEs(:,:,:,:,:,1:end));
end

end
