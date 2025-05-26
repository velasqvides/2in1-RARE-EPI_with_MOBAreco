% set the paths for BART and the functions required for image reconstruction
if isunix
bartPath = 'tools/bart_v09/bart';
run(fullfile(bartPath, 'startup.m'));
end

codePath = fullfile(pwd, 'code');
addpath(genpath(codePath));

toolsPath = fullfile(pwd, 'open_source_tools');
addpath(genpath(toolsPath));

previousDir = pwd;
cd('open_source_tools/arrShow');
arrShow.registerPaths();
cd(previousDir);

clear codePath bartPath previousDir;