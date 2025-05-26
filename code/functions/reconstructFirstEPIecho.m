function reconstructFirstEPIecho(kSpaceEPI,trajEPI,protPara)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here
% kSpacefirstEchoEPI = zeros(1,protPara.baseRes*protPara.oversamplingFactor,protPara.nSpokes,protPara.nChannels);
% tmp = kSpaceEPI(:,:,:,:,1,1,1);
% kSpacefirstEchoEPI(1,:,:,:) = tmp(1,:,:,:);
% trajFirstEchoEPI = trajEPI(:,:,:,1,1,1);
% trajFirstEchoEPI = squeeze(trajFirstEchoEPI);
image = RecontructImageGridding(kSpaceEPI,trajEPI,protPara);
as(image)
end