
function [kSpace, gdelcalib, protPara, config] = extractRecoInfo(measFileLocation,removeOversampling)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%% 0. initialize a struct variable called config
[filepath,name,ext] = fileparts(measFileLocation);
config.filename = strcat(name,ext);
config.pathname = filepath;
config.preprocessing.REMOVE_OVERSAMPLING = removeOversampling;
processSiemensRaw;
% calibrationScans = gdelcalib;
% kSpaceData = kSpace;
%% 1. Extract required recoonstruction parameters
wipPara = extractWipPara(MrProt); % extract parameters from special card
protPara.ETL_RARE = wipPara.lNEchoesRARE;
protPara.ETL_EPI = wipPara.lNEchoesEPI;
protPara.ETLfull = protPara.ETL_RARE + protPara.ETL_EPI;
protPara.interleaved = wipPara.interleaved;
protPara.bSameSpokeInET = wipPara.bSameSpokeInET;
protPara.ViewOrdering = wipPara.ViewOrdering;
protPara.ESP_RARE = wipPara.dESP_RARE;
protPara.ESP_EPIfirst = wipPara.dESP_EPIFirst;
protPara.ESP_EPI = wipPara.dESP_EPI;
if ~isempty(wipPara.lNdummies)
    protPara.nDummies = wipPara.lNdummies;
else
    protPara.nDummies = 0;
end
protPara.nSpokesRAREall = MrProt.sKSpace.lRadialViews * protPara.ETL_RARE;
protPara.nSpokesEPIall = MrProt.sKSpace.lRadialViews * protPara.ETL_EPI;
protPara.baseRes = MrProt.sKSpace.lBaseResolution; 
protPara.nSpokes = MrProt.sKSpace.lRadialViews;
protPara.nSpokesFull = protPara.nSpokes * protPara.ETLfull;
nChannels = size(MrProt.sCoilSelectMeas.aRxCoilSelectData{1}.asList,2);
protPara.nChannels = nChannels;
protPara.nSlices = MrProt.sSliceArray.lSize;
oversamplingFactor = size(kSpace,1) / protPara.baseRes;
protPara.oversamplingFactor = oversamplingFactor;
protPara.ucTrajectory = MrProt.sKSpace.ucTrajectory;


%% 3. Get data dimensions in fixed order [Col Lin Cha Sli]
if (exist('kSpace','var') && protPara.nSlices > 1)
    kSpace = permute(kSpace,[1 2 4 3]);
end % end if

if ~(exist('gdelcalib','var') )
    gdelcalib = [];
end % end if

%% 4. Save required variables
config.pathname = append(filepath,"/");
folders = split(config.pathname,"/");
subjectName = folders{end-1,1};
% timeStamp = strrep(strrep(datestr(now),' ','_'),':','-');
% subjectName = append(subjectName,'_',timeStamp);
measurementName = config.filename;
[~,measurementName,ext] = fileparts(measurementName);
folderName = '~/Documents/projects/MBR_automatized/reconstructions';
dirToSave = fullfile(folderName,subjectName,measurementName);
mkdir(dirToSave);
config.dirToSave = dirToSave;
end % end function