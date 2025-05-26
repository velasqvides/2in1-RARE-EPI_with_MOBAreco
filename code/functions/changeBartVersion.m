function changeBartVersion(version)
if version == 7
    functionDir = fileparts(mfilename('fullpath'));
    bartPath = fullfile(functionDir,'../../tools/bart_v07/bart');
    run(fullfile(bartPath, 'startup.m'));
elseif version == 9
    functionDir = fileparts(mfilename('fullpath'));
    bartPath = fullfile(functionDir,'../../tools/bart_v09/bart');
    run(fullfile(bartPath, 'startup.m'));
end

end


