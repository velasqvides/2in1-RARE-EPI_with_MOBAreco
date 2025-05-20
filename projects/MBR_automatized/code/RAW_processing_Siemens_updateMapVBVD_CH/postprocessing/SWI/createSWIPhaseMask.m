function [ phaseMask ] = createSWIPhaseMask( phaseImage, maskType )
%CREATESWIPHASEMASK creates a linear SWI phase mask
%   Usage:
%   [ phaseMask ] = createSWIPhaseMask( phaseImage, maskType )
%   
%   maskType    = expects a string 'negative' or 'positive'
%   phaseImage  = phase image with phase in the max range of -pi to +pi
%

%   2015-03-05, Till Huelnhagen, greatly improved speed and memory
%               consumption by avoiding loop

%% create phase mask for SWI
phaseMask=ones(size(phaseImage));
if strcmp(maskType,'positive')
    phaseMask(phaseImage > 0) = 1 - phaseImage(phaseImage > 0) ./ pi;
else
    phaseMask(phaseImage < 0) = 1 + phaseImage(phaseImage < 0) ./ pi;
end

end

