function [wipPara] = extractWipPara(MrProt)
%This function extracts parameters stored in MrProt.WipMemBlock
%   INPUT:   MrProt or MeasYaps struct
%   OUTUPT:  wipPara struct with fieldnames according to parameter names
%
%   CAUTION: If parameters are added/removed from special card, names need
%            to be updated accordingly

wipPara = [MrProt.sWipMemBlock.adFree MrProt.sWipMemBlock.alFree];

wipParaNames = {'dRFbwtpExc','dRFbwtpRef','dRFDurExc','dRFDurRef','dESP_RARE','dESP_EPIFirst','dESP_EPI',...
                 'lNEchoesRARE', 'lNEchoesEPI', 'lNGdcScans',...
                 'bSliCrush', 'bGSpoilDirVarRO', 'tinyGA_RARE','tinyGA_EPI', 'ViewOrdering'}; 

if (size(wipParaNames,2) == size(wipPara,2))
    wipPara = cell2struct(wipPara,wipParaNames,2);
else
    wipPara = [];
    fprintf('\n WARNING: Number of names does not match number of parameters in WipMemBlock. Skipping extraction of WipMemBlock... \n Check for changes in special card or outdated raw data. Probably adjust naming manually. \n')
end

end

