function [T2, R2, M0, sensRARE, synthesizedRAREimages, binaryMaskRARE_] = prepareRAREdata(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages, protPara)
T2 = prepareMaps(T2, binaryMaskRARE, protPara);
R2 = prepareMaps(R2, binaryMaskRARE, protPara);
M0 = prepareGrayImages(M0, binaryMaskRARE, protPara);
sensRARE = prepareImages(sensRARE, binaryMaskRARE, protPara);
synthesizedRAREimages = prepareImages(synthesizedRAREimages, binaryMaskRARE, protPara);
binaryMaskRARE_ = prepareGrayImages(binaryMaskRARE, binaryMaskRARE, protPara);
end

