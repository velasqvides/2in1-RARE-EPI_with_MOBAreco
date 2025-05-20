function [ T2starMap ] = T2StarMappingARLO( combinedImage,MrProt, newDataDims )
%T2STARMAPPINGARLO Summary of this function goes here
%   Detailed explanation goes here
% 2014-08-05, created by Till Huelnhagen

% delete channel dimension if channels have been combined
newDataDims(find(strcmp(newDataDims,'Cha')))=[];

% verify multiecho scan
ecoDim = find(strcmp(newDataDims,'Eco'));
if isempty(ecoDim)
    disp('ERROR: No multiecho scan. Multiecho data required for T2* mapping');
    return
end

% verify minimum echo number of 3
if size(combinedImage,ecoDim) < 2
    disp('ERROR: Not enough echoes for ARLO T2* mapping. Minimum of three echoes requiered');
    return
end

% get magnitude
img = abs(combinedImage);

% make echo dimension last in array
if ecoDim < numel(newDataDims)
    shiftOrder=[1:ecoDim-1,ecoDim+1:numel(newDataDims), ecoDim];
    img=permute(img,shiftOrder);
end

% calculate deltaTE in seconds
TE=(MrProt.TE/1E6)';
TE=TE(1:MrProt.Contrasts);
deltaTE=TE(2)-TE(1);
% TEdiff=circshift(TE,[1,0]);
% TEdiff(1)=0;

T2starMap=ARLOT2starFit(img, deltaTE);

fprintf('\nARLO T2* mapping finished\n');

end

