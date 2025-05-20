function kSpace_cc = performCoilCompression(kSpace,nVirtualCoils)
%% 9 coil compression
% cc does nto work in BART v0.9.00 so i have to go back to v0.7.00 for cc
% functionDir = fileparts(mfilename('fullpath'));
% bartPath = fullfile(functionDir,'../../../../tools/bart_v07/bart');
% run(fullfile(bartPath, 'startup.m'));
changeBartVersion(7)
kSpace_cc = bart(sprintf('cc -p%i -A -S',nVirtualCoils),kSpace);
changeBartVersion(9)
% return to bart v0.9.00
% bartPath = fullfile(functionDir,'..filesep..filesep..filesep..fileseptools/bart_v09/bart');
% run(fullfile(bartPath, 'startup.m'));
end