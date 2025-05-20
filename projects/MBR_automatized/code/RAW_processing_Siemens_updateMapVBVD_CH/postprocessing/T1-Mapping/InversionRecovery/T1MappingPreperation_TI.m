%% T1 Mapping Preperation
% Thomas Eigentler
% 18-10-2018


%% Settings

clear all

addpath(genpath('K:\MATLAB\SVN_repository\trunk\RAW_processing')); % Add non linear FFT Code to Matlab search pat

pathname = 'P:\Measurement_7T\20181020_T1T2_IkeaPhantom';

%% Add Images

index = 1;

while 1
    clear config
    config.pathname = pathname;
    processSiemensRaw
    
    Temp.combinedImage(:,:,index) = combinedImage(:,:);
    Temp.TI(index) = MrProt.TI;
    
    reply = input('To add a File select 1, to stop press [Enter]');
    % prevent window from disappearing
    
    index = index + 1;
    %in case [Enter] is pressed
    reply(isempty(reply))=0;
    if (reply > 0  || reply < 0)
        continue
    else
        break
    end  
end

MrProt.TI = Temp.TI;
combinedImage = Temp.combinedImage;

