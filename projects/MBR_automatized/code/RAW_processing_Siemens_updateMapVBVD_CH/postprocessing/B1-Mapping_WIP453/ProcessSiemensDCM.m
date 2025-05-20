function [dDcm, dDcmInfo, filenames] = ProcessSiemensDCM(config)

olddir = cd;

if not(isfield(config,'filename')) && isfield(config,'pathname')
    [config.filename, config.pathname, ~] = uigetfile( ...
        {'*.dcm','Meas files VB17 (*.)'}, ...
        'Pick a Meas file', config.pathname);
elseif not(isfield(config,'filename'))
    [config.filename, config.pathname, ~] = uigetfile( ...
        {'*.dcm','Meas files VB17 (*.dcm)'}, ...
        'Pick a Meas file', 'meas.dcm');
end

if isempty(config.filename), disp('ERROR: No Measurement data are found.'); return; end

filenameDcm = fullfile(config.pathname, config.filename);

dDcm(:,:) = dicomread(char(filenameDcm));
dDcmInfo = dicominfo(char(filenameDcm)); 
filenames = filenameDcm;
dSliceLocationTmp = dDcmInfo.SliceLocation;
  
cd(olddir);
end

