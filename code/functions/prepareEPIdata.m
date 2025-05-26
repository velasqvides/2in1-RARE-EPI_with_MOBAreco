function [T2star, R2star, M0star, sensEPI, synthesizedEPIimages, B0, initB0, initMaps] = prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE, sensEPI, synthesizedEPIimages, B0, initMaps, protPara)
T2star = prepareMaps(T2star, binaryMaskRARE, protPara);
R2star = prepareMaps(R2star, binaryMaskRARE, protPara);
M0star = prepareGrayImages(M0star, binaryMaskRARE, protPara);
sensEPI = prepareImages(sensEPI, binaryMaskRARE, protPara);
synthesizedEPIimages = prepareImages(synthesizedEPIimages, binaryMaskRARE, protPara);
B0 = prepareMaps(B0, binaryMaskRARE, protPara);
initMaps = prepareInitMaps(initMaps, binaryMaskRARE, protPara);
initB0 = initMaps(:,:,4,:);
initB0 = squeeze(initB0 );
end