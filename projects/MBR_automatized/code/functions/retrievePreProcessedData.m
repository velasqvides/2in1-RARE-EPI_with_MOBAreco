function [kSpaceRARE, kSpaceEPI, trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = retrievePreProcessedData(folderWithPreProcessedData)
cd(folderWithPreProcessedData);
kSpaceRARE = readcfl('kSpaceRARE');
kSpaceEPI  = readcfl('kSpaceEPI');
trajRARE   = readcfl('trajRARE');
trajEPI    = readcfl('trajEPI');
TEsRARE    = readcfl('TEsRARE');
TEsEPI     = readcfl('TEsEPI');
load config.mat config;
load protPara.mat protPara;
end