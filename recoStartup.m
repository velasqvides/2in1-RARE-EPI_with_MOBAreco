% set the paths for BART and the functions required for image reconstruction
mainFolder = fileparts(mfilename('fullpath'));
bartPath = fullfile(mainFolder, 'open_source_tools', 'bart_v09', 'bart');
run(fullfile(bartPath, 'startup.m'));

codePath = fullfile(mainFolder, 'code');
addpath(genpath(codePath));

mapVBVDPath = fullfile(mainFolder, 'open_source_tools', 'RAW_processing_Siemens');
addpath(genpath(mapVBVDPath));

colorResourcesPath = fullfile(mainFolder, 'open_source_tools', 'mfuderer-colorResources-37639c3');
addpath(genpath(colorResourcesPath));


clearvars;

