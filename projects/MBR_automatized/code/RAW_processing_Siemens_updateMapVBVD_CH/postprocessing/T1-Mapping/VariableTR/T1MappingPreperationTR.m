%% T1 Mapping Preperation
% Thomas Eigentler
% 18-10-2018


%% Settings

%clear all;

%addpath(genpath('K:\MATLAB\SVN_repository\trunk\RAW_processing')); % Add non linear FFT Code to Matlab search pat

%pathname = 'P:\Measurement_7T\20181020_T1T2_IkeaPhantom';

%% Add Images


pathMainDir = uigetdir('O:\Users\Carl\Radial_RARE-EPI\rawData\');
filenames = dir([pathMainDir '\*.dat']);
config.pathname = pathMainDir;
for j=1:length(filenames)
    % Carl Herrmann, 19-10-2019
    % remove specification of config here, instead config should be adjusted in
    % individual matlab scripts, rather than a general code repository
    %config.pathname = pathname;
    config.filename = filenames(j).name;
    processSiemensRaw;
    Temp.combinedImage(:,:,j) = combinedImage(:,:);
    % addapt to Software version VE11 (i.e. TR -> alTR (cell array) )
    Temp.TR(j) = cell2mat(MrProt.alTR);
    
%     reply = input('To add a File select 1, to stop press [Enter]');
    % prevent window from disappearing
       
%     %in case [Enter] is pressed
%     reply(isempty(reply))=0;
%     if (reply > 0  || reply < 0)
%         continue
%     else
%         break
%     end  
end

MrProt.alTR = Temp.TR;
combinedImage = Temp.combinedImage;
