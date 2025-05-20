%% PRFS main funktion
% Thomas Eigentler
% 2019-08-05

%% Revision Log

% 2019-08-05: Thomas Eigentler: Generation

%% Initialize

addpath(genpath('K:\MATLAB\SVN_repository\trunk\RAW_processing_Siemens'));

%% MRThermometry

DataFolder = 'P:\_Project-ThermalMR_Haopeng_2Channel_SignalGen\7TMeasurement';

[FileNamePre, PathnamePre] = uigetfile('*.dat','pick scan data prior heating',DataFolder);
imgPre = [PathnamePre, FileNamePre];
[FileNamePost, PathnamePost] = uigetfile('*.dat','pick scan data post heating',PathnamePre);
imgPost = [PathnamePost, FileNamePost];

tempMaps = calcTempMaps(imgPre,imgPost,'driftCorrectionMode','mask','smoothKer',1);
%tempMaps =
%calcTempMaps(imgPre,imgPost,'driftCorrectionMode','mask','smoothKer',ones(3,1));
%with smoothing!!!

close all

%% Visualization And Data Save

% Generate Image and Save it
MaxTemp = 5;
FontSize = 20;
FigureSize = 800;

WriteImage = true;
ExportPath = [PathnamePre, 'Processed_MRThermomety\'];
mkdir(ExportPath);

for Slice = 1:size(tempMaps.TemperatureMap,3)
    for Image = 1:size(tempMaps.TemperatureMap,4)
        ExportFileName = ['TemperatureMap-Slice', num2str(Slice), '-Image', num2str(Image)];
        TemperaturePlotter(tempMaps.TemperatureMap(:,:,Slice,Image), ExportPath, ExportFileName, MaxTemp, FontSize, FigureSize)
        close all
    end
end

% Save Data
exp = [ExportPath '\TemperatureMap.mat'];
save(exp,'tempMaps','-v7.3');

