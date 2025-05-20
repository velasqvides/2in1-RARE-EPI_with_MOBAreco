function [image, visuPars] = importBrukerReco(filepath)
% pvtools must be on the searchpath

Visu = readBrukerParamFile([filepath 'visu_pars']);

[image, visuPars] = readBruker2dseq([filepath 'pdata' filesep '1' filesep '2dseq'], Visu);

end
