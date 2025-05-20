% set the paths for BART and the functions required for image reconstruction
if isunix
bartPath = '../../tools/bart_v09/bart';
run(fullfile(bartPath, 'startup.m'));
end
codePath = fullfile(pwd, 'code');
addpath(genpath(codePath));

previousDir = pwd;
cd('../../tools/arrShow');
arrShow.registerPaths();
cd(previousDir);

clear codePath bartPath previousDir;