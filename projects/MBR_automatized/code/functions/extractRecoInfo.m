function [kSpace, gdelcalib, protPara, config] = extractRecoInfo(measFileLocation, removeOversampling, applyNoiseDecorrelation)
%% 0. initialize a struct variable called config
inputPath = char(measFileLocation);
[filepath,name,ext] = fileparts(inputPath);
config.filename = strcat(name,ext);
config.pathname = filepath;
config.preprocessing.REMOVE_OVERSAMPLING = removeOversampling;
config.preprocessing.IGNORE_SEGMENTS = 1;
config.preprocessing.APPLY_NOISE_DECORRELATION = applyNoiseDecorrelation;
config.preprocessing.READ_ADDITIONAL_SCANS = 1;
config.preprocessing.DO_AVERAGE = 0;
config.preprocessing.DO_SVDBasedChCombination = 0;
processSiemensRaw;
if (MrProt.sKSpace.ucTrajectory == 2) % 2 -> radial
scaleFactor = 10000;
integerPartOfAngles = icePara(1, :);
decimalPartOfAngles = icePara(2, :) / scaleFactor;
allAngles = integerPartOfAngles + decimalPartOfAngles;
protPara.allAngles = allAngles;
end
% calibrationScans = gdelcalib;
% kSpaceData = kSpace;
%% 1. Extract required recoonstruction parameters

protPara.baseRes = MrProt.sKSpace.lBaseResolution; 
oversamplingFactor = size(kSpace,1) / protPara.baseRes;
protPara.oversamplingFactor = oversamplingFactor;
nChannels = size(MrProt.sCoilSelectMeas.aRxCoilSelectData{1}.asList,2);
protPara.nChannels = nChannels;
protPara.nSlices = MrProt.sSliceArray.lSize;
protPara.trajectoryType = MrProt.sKSpace.ucTrajectory;
protPara.sliceThickness = MrProt.sSliceArray.asSlice{1, 1}.dThickness;
if (MrProt.sKSpace.ucTrajectory == 2) % 2 -> radial
protPara.bSameSpokeInET=1; %required for old gradient correcton function
wipPara = extractWipPara(MrProt); % extract parameters from special card
protPara.ETL_RARE = wipPara.lNEchoesRARE;
protPara.ETL_EPI = wipPara.lNEchoesEPI;
protPara.ETLfull = wipPara.lNEchoesRARE + wipPara.lNEchoesEPI;

protPara.ViewOrdering = wipPara.ViewOrdering;
protPara.ESP_RARE = wipPara.dESP_RARE;
protPara.ESP_EPIfirst = wipPara.dESP_EPIFirst;
protPara.ESP_EPI = wipPara.dESP_EPI;
protPara.nSpokesRAREall = MrProt.sKSpace.lRadialViews * protPara.ETL_RARE;
protPara.nSpokesEPIall = MrProt.sKSpace.lRadialViews * protPara.ETL_EPI;
protPara.nSpokes = MrProt.sKSpace.lRadialViews;
protPara.nSpokesFull = protPara.nSpokes * protPara.ETLfull;
protPara.tinyGA_RARE = wipPara.tinyGA_RARE;
protPara.tinyGA_EPI = wipPara.tinyGA_EPI;
else % 1 -> cartesian
protPara.PhaseEncodingLines = MrProt.sKSpace.lPhaseEncodingLines;
protPara.ETL = MrProt.lContrasts;
protPara.ESP = MrProt.alTE{1}/1000;
end
%% 3. Get data dimensions in fixed order [Col Lin Cha Sli]
if (exist('kSpace','var') && protPara.nSlices > 1)
    kSpace = permute(kSpace,[1 2 4 3 5]);
end % end if

if ~(exist('gdelcalib','var') )
    gdelcalib = [];
end % end if

%% 4. Save required variables
config.pathname  = append(filepath,filesep);
folders          = split(config.pathname,filesep);
subjectName      = folders{end-1,1};
functionDir      = fileparts(mfilename('fullpath'));
folderName       = fullfile(functionDir, '..', filesep, '..', filesep, 'reconstructions');       
dirToSave        = fullfile(folderName, subjectName,name);
mkdir(dirToSave);
config.dirToSave = dirToSave;
end % end function