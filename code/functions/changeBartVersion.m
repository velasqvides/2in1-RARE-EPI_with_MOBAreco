function changeBartVersion(version)

functionDir = fileparts(mfilename('fullpath'));

if version == 7
    bartPath = fullfile(functionDir,'../../open_source_tools/bart_v07/bart');
    run(fullfile(bartPath, 'startup.m'));
elseif version == 9
    bartPath = fullfile(functionDir,'../../open_source_tools/bart_v09/bart');
    run(fullfile(bartPath, 'startup.m'));
else
    error('Unsupported BART version: %d. Supported versions are 7 and 9.', version);
end

end


