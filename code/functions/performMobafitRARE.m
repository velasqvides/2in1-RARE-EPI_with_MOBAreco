function [M0, R2, T2] = performMobafitRARE(img,TEs,iter)
nSlices = size(img,7);
nSamples1 = size(img,1);
nSamples2 = size(img,2);
mobaFitMaps = zeros(nSamples1,nSamples2,2,nSlices);
for i=1:nSlices
maps = bart( sprintf( 'mobafit -T -i%i',iter ), TEs, img(:,:,:,:,:,:,i) );
maps = squeeze(maps);
mobaFitMaps(:,:,:,i) = maps;
end
M0 = squeeze(mobaFitMaps(:,:,1,:));
R2 = squeeze(mobaFitMaps(:,:,2,:));
T2 = 1./R2;
T2 = T2.*1000;
end