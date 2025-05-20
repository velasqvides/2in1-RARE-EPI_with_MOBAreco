function changeBartVersion(version)

if ispc
    if version == 7
        toolbox_path = '${HOME}/bart_v07/bart/bart';
    elseif version == 9
        toolbox_path = '${HOME}/bart_v09/bart/bart';
    end
    cmd1 = 'wsl rm /usr/local/bin/bart';
    [~,~] = system(cmd1);
    cmd4 = ['wsl ln -s ' toolbox_path ' /usr/local/bin/bart'];
    [~,~] = system(cmd4);
elseif isunix
    if version == 7
        functionDir = fileparts(mfilename('fullpath'));
        bartPath = fullfile(functionDir,'../../../../tools/bart_v07/bart');
        run(fullfile(bartPath, 'startup.m'));
    elseif version == 9
        functionDir = fileparts(mfilename('fullpath'));
        bartPath = fullfile(functionDir,'../../../../tools/bart_v09/bart');
        run(fullfile(bartPath, 'startup.m'));
    end

end

end
