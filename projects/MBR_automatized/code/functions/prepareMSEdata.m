function [M0MSE, R2MSE, T2MSE, imagesMSE] = prepareMSEdata(M0MSE, R2MSE, T2MSE, imagesMSE,  binaryMaskRARE, protPara)
T2MSE = prepareMaps(T2MSE, binaryMaskRARE, protPara);
R2MSE = prepareMaps(R2MSE, binaryMaskRARE, protPara);
M0MSE = prepareGrayImages(M0MSE, binaryMaskRARE, protPara);
imagesMSE = prepareImages(imagesMSE, binaryMaskRARE, protPara);
end