function [M0star, T2star, R2star, B0] = performMobafitEPI(img,TEs,iter)
nSlices = size(img,7);
nSamples1 = size(img,1);
nSamples2 = size(img,2);
mobaFitMaps = zeros(nSamples1,nSamples2,3,nSlices);
for i=1:nSlices
maps = bart( sprintf( 'mobafit -G -m3 -i%i',iter ), TEs, img(:,:,:,:,:,:,i) );
maps = squeeze(maps);
mobaFitMaps(:,:,:,i) = maps;
end
% as(mobaFitMapsRARE)
M0star = mobaFitMaps(:,:,1,:);
B0 = mobaFitMaps(:,:,3,:);
R2star = mobaFitMaps(:,:,2,:);
T2star = 1./R2star;
T2star = T2star.*1000;
T2star = squeeze(T2star);
M0star = squeeze(M0star);
R2star = squeeze(R2star);
B0 = squeeze(B0);
end




