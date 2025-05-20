function TempMap  = calcTempMaps(inputFile_preheat, inputFile_postheat, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% parse input
allowedDriftCorrModes = {'none','mask','val'};
p = inputParser;
addParameter(p,'b0FieldStrength',6.98,@(x) isscalar(x) && x>0);
addParameter(p,'nucleusGamma',267.5221900e6,@(x) isscalar(x));
addParameter(p,'tempAlpha',-0.0106*1.0E-6,@(x) isscalar(x));
addParameter(p,'smoothKer',gausswin(5,1),@(x) isvector(x));
addParameter(p,'maskThresh',0.05,@(x) isscalar(x));
addParameter(p,'driftCorrectionMode','none', @(x) any(contains(allowedDriftCorrModes, x)));
addParameter(p,'driftCorrectionValue',0, @(x) isscalar(x));

parse(p, varargin{:});

b0Strength = p.Results.b0FieldStrength;
smoothKer = p.Results.smoothKer;
maskThresh = p.Results.maskThresh;
driftCorrectionMode = p.Results.driftCorrectionMode;
driftCorrectionValue = p.Results.driftCorrectionValue;
nucleusGamma = p.Results.nucleusGamma;
tempAlpha = p.Results.tempAlpha;

% generate config for processSiemensRaw
config.preprocessing.REMOVE_OVERSAMPLING = 1;
[inputFile_preheat] = GetFullPath(inputFile_preheat);
[config.pathname, inputFile_preheat, inputExt] = fileparts(inputFile_preheat);
config.filename = [inputFile_preheat inputExt];

% process raw data for preheating
processSiemensRaw;
combinedImage = double(combinedImage);

% get resolutions
TempMap.dx = MrProt.SliceArray.Slice(1).ReadoutFOV/size(combinedImage,1);
TempMap.dy = MrProt.SliceArray.Slice(1).PhaseFOV/size(combinedImage,2);

if ndims(combinedImage) == 3
    combinedImage_preheat = permute(combinedImage,[1 2 4 3]); % make 2D into 3D with 1 slice
else
    combinedImage_preheat = combinedImage;
end
MrProt_preheat = MrProt;

% generate config for processSiemensRaw
clear config;
config.preprocessing.REMOVE_OVERSAMPLING = 1;
[inputFile_postheat] = GetFullPath(inputFile_postheat);
[config.pathname, inputFile_postheat, inputExt] = fileparts(inputFile_postheat);
config.filename = [inputFile_postheat inputExt];

% process raw data for postheating
processSiemensRaw;
combinedImage = double(combinedImage);

if ndims(combinedImage) == 3
    combinedImage_postheat = permute(combinedImage,[1 2 4 3]); % make 2D into 3D with 1 slice
else
    combinedImage_postheat = combinedImage;
end
MrProt_postheat = MrProt;

clear sensitivity kSpace combinedImage

TempMap.imagePre = combinedImage_preheat;
TempMap.imagePost = combinedImage_postheat;

% smoothing
ker = smoothKer;
ker = ker*ker';
ker = ker./sum(ker(:));
combinedImage_preheat = convn(combinedImage_preheat,ker,'same');
combinedImage_postheat = convn(combinedImage_postheat,ker,'same');

TempMap.imagePreSmoothed = combinedImage_preheat;
TempMap.imagePostSmoothed = combinedImage_postheat;

% masking calc
maskImage = abs(combinedImage_preheat(:,:,:,1));
maxVal = max(maskImage(:));
mask = maskImage >= maskThresh*maxVal;

% complex phase differences for all echos
dPhi = angle(combinedImage_postheat./combinedImage_preheat);
dPhi = dPhi(:,:,:,2:end); % discard first echo

% drift correction
switch driftCorrectionMode
   
    case 'mask'
                imt = imtool3D(abs(TempMap.imagePre(:,:,:,1)));
                input('Press return to continue...');
                TempMap.driftMask = imt.getMask;
                tmpDriftImage = dPhi(:,:,:,1); %use first echo for drift correction (lowest TE - best SNR)
                TempMap.driftVal = mean(tmpDriftImage(TempMap.driftMask));
                
    case 'val'
                TempMap.driftVal = driftCorrectionValue;
                
    case 'none'
                TempMap.driftVal = 0;
                
end

% drift correction & use mask
dPhi = (dPhi - TempMap.driftVal);
TEDiff = (MrProt_preheat.TE(2:MrProt.Contrasts)'-MrProt_preheat.TE(1))*1e-6;

TempMap.TemperatureMap = dPhi./(b0Strength*nucleusGamma*tempAlpha*permute(TEDiff,[4 3 2 1]));
TempMap.TemperatureMapMasked = TempMap.TemperatureMap.*mask;
TempMap.Mask = mask;

end