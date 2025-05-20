function [B0Map_Hz,detailedB0postprocessingData]=B0Mapping_20140923(combinedImage,MrProt,oldDataDims,config)
% This function calculates B0map of a multi-GRE-sequence
% function uses the B1Mapping script as blueprint
%
% INPUT:
% combinedImage = complex images with at least 2 gradient echoes
% MrProt = Struct with MR protocol parameters
% oldDataDims =  cell with description of dimensions of kSpace data
% config = struct with information such as filename, pathname and savingString
%          if no saving string is defined, no data will be saved
%          if config.postprocessing.B0_MAPPING_savingSTRING = 'someString' all postprocessing data will be saved in someString_MIDxx.mat
%          format of config is determined by processSiemensRaw script!
%
% OUTPUT:
% B0Map_Hz = B0 map in [Hz]
% detailedB0postprocessingData = structure of other relevant data
%
% modifications by Olli Kraus
% 2014.10.01

%% first we need to transform the data to double, since all other function
% use doubles. We would need to rewrite the whole program, especially the
% sensitivity calculation, otherwise
combinedImage=cast(combinedImage,'double');

Info.usedMfile=[mfilename('fullpath'),'.m'];
Info.TimeOfExecution=datestr(clock);

Info.usedSeq=MrProt.SequenceFileName;

Info.RawPath=config.pathname;
Info.MID=regexpi(config.filename,'_','split');
Info.MID=Info.MID{2};

Info.ProcessedFolder=strcat(Info.MID ,'_Processed_',datestr(now,'yyyymmdd'));
Info.ProcessedPath=fullfile(config.pathname,Info.ProcessedFolder);

if isempty(config.postprocessing.B0_MAPPING_savingSTRING)
    Info.SavingString=0;
else
    Info.SavingString=genvarname(config.postprocessing.B0_MAPPING_savingSTRING);
    mkdir(Info.ProcessedPath);
end

%% delete 'Channels' from the dimensions because the channels were already combined
oldDataDims(strcmp(oldDataDims,'Cha'))=[];
%% permute the Eco-Dimension on position 3
EcoPos=find(~cellfun(@isempty,strfind(oldDataDims,'Eco')));
oldOrder=1:length(oldDataDims);
oldOrderStart=oldOrder(1:2);
oldOrderMid=oldOrder(3:EcoPos-1);
oldOrderEnd=oldOrder(EcoPos+1:end);

permOrder=[oldOrderStart,EcoPos,oldOrderMid,oldOrderEnd];
% for i=1:length(oldDataDims)
%     Info.DataDims{i}=oldDataDims{permOrder(i)};
% end
clear oldOrder*
%permute the data according to the new order
CmplxImage=permute(combinedImage,permOrder);
MagMax=max(abs(CmplxImage(:)));
CmplxImage=CmplxImage/MagMax; % normalize magnitude to 1
% extract dimensions into variables
nPar=size(CmplxImage,4); % number of 3D partitions/slices
% number of echoes
Info.noEchoes=size(CmplxImage,3);
if Info.noEchoes<2
    disp('Error: no multiple echoes detected -> end of script! ')
    return
else
    fprintf('%d echoes detected\n',Info.noEchoes);
end

%% %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%           %%%
%%%  MASKING  %%%
%%%           %%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
Magnitude1=squeeze(abs(CmplxImage(:,:,1,:)));
% this method is robust but can't cut away inner parts (eg lung tissue)
%[Masks.Mask,Masks.Perim,Masks.Image]=masking_best_EM3D(Magnitude1);

%this method can also cut away parts within an object
[Masks.MagThres,Masks.Image,Info.relativeMagnitudeThreshold]=masking_by_thresholding_3D(Magnitude1);
Masks.MagThres=logical(Masks.MagThres);
close(gcf)
for cPar=1:nPar
    Image_plot3D(1,nPar,cPar,0,1,'jet',Masks.MagThres(:,:,cPar),0,1)
end
set(gcf, 'name', 'ThresholdMask')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'ThresholdMask'),'-dtiff')
end
close(gcf)
%% optional: draw ROI for every slice
Info.enROI=input('Do you want to draw a ROI [1] or not [Enter]? ');
%Info.enROI=1;
%close all
if  (~isempty(Info.enROI))
    for cPar=1:nPar
        figure(11)
        imagesc(Masks.Image(:,:,cPar));
        axis image
        figID=gca;
        ROI=imfreehand(figID);
        Masks.lineOfROI{cPar}=getPosition(ROI);
        %ROI=imellipse(figID);
        %ROI=imrect(figID);
        %wait(ROI);
        Masks.ROI(:,:,cPar)=createMask(ROI);
        Image_plot3D(16,nPar,cPar,0,1,'jet',Masks.ROI(:,:,cPar).* Masks.Image(:,:,cPar),0,1)
    end
    set(gcf, 'name', 'ROIMask')
    if Info.SavingString
        Image_save(fullfile(Info.ProcessedPath,'ROIMask'),'-dtiff')
    end
else
    Info.enROI=0;
    disp('no ROI drawn')
end
close all

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  B0 CALCULATIONS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate phase differences for B0 map
% "phasedifference_B0 = phase_Echo2 - phase_Echo1 = phase_Echo4 - phase_Echo3"

%imPhaDiffB0=squeeze(mod(angle(CmplxImage(:,:,2,:))-angle(CmplxImage(:,:,1,:))+pi,2*pi)-pi);
imPhaDiffB0=squeeze(phaseDiff(CmplxImage(:,:,2,:),CmplxImage(:,:,1,:)));

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%%  PHASE UNWRAPPING  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nPar~=1 % multislice or 3D Maps
    for cPar=1:nPar % display B0 phase differences to check for phase wraps
        data2plot=imPhaDiffB0(:,:,cPar).*Masks.MagThres(:,:,cPar);
        Min=min(data2plot(:));
        Max=max(data2plot(:))+0.01;
        Image_plot3D(21,nPar,cPar,Min,Max,'jet',data2plot,0,1)
    end
    disp('Please check B0 phase differences for phase wraps!');
    Info.enUnwrapping=input('Does your data need phase unwrapping [1] or not [Enter]? ');
else
    Info.enUnwrapping=1;
end
close all
if (~isempty(Info.enUnwrapping))
    disp('===>> Unwrapping, please wait...')
    for cPar=1:nPar
             inpt=1;   
%        imagescwithnan(imPhaDiffB0(:,:,cPar)./Masks.MagThres(:,:,cPar),jet(256),[1 1 1])
%        set(gca,'Title',text('String',['Phase difference. slice: ', num2str(cPar)]));
%        axis image
%        colormap(jet)
%        inpt=input('Are there phase wraps [any number] or none [Enter]? ');
        if (~isempty(inpt))
            cmplx1=Masks.MagThres(:,:,cPar).*Magnitude1(:,:,cPar).*exp(1i.*imPhaDiffB0(:,:,cPar));
            temp=angle(cmplx1);
            if max(temp(:))~=min(temp(:))
                while 1
                    close all
                    imagesc(temp)
                    colorbar
                    axis image
                    set(gca,'Title',text('String',['slice: ',num2str(cPar),', Pick a point, where \Delta B_0 is small and with good signal!']))
                    [xpoint,ypoint] = ginput(1);
                    % perform unwrapping
                    imPhaDiffB0(:,:,cPar)=Unwrap2D_QualityGuided_inputXY(cmplx1, Info.relativeMagnitudeThreshold ,0,xpoint,ypoint);
                    imagesc(imPhaDiffB0(:,:,cPar))
                    colormap jet
                    set(gca,'Title',text('String',['Phase difference used for B_0-Map, slice: ' num2str(cPar)]))
                    axis image
                    %colorbar
                    
                    inpt2=[];
                    %inpt2=input('DONE! Happy with the result [any key] or not [1]? ');
                    inpt2(isempty(inpt2))=0;
                    if inpt2==0
                        break
                    end
                end
            end
        end
    end
    Info.UnwrappingInput=inpt;
    disp('===>> Unwrapping done!')
end


%% calculate B0 Map
disp('===>> Calculating B0 maps...');
B0.DeltaTE12_ms=(MrProt.TE(2)-MrProt.TE(1))*1e-3;
disp(['DeltaTE12_ms = ',num2str(B0.DeltaTE12_ms),' ms']);

B0.Map_Hz=imPhaDiffB0/2/pi/ B0.DeltaTE12_ms .*Masks.MagThres*1e3;
Info.MaxDeltaB0=max(abs(B0.Map_Hz(:)));
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-Info.MaxDeltaB0,Info.MaxDeltaB0,jet,B0.Map_Hz(:,:,cPar)./Masks.MagThres(:,:,cPar),1,0)
end
set(gcf, 'name', 'B0 map')

%% cut off outliers
disp(['Maximum B0 is ',num2str(round(Info.MaxDeltaB0)),' Hz']);
Info.B0cut=input('choose cut off value ([Enter] to use max)> ');
if isempty(Info.B0cut)
    Info.B0cut=Info.MaxDeltaB0;
end
B0.B0cut=Info.B0cut;
B0.outliers=logical(abs(B0.Map_Hz)>Info.B0cut);
B0.noOutliers=sum(B0.outliers(:));
Masks.B0outliers=B0.outliers;

%% save image and data
for cPar=1:nPar
    Image_plot3D(3,nPar,cPar,-Info.B0cut,Info.B0cut,'jet',B0.Map_Hz(:,:,cPar).*(1-B0.outliers(:,:,cPar)),1,1)
end
set(gcf, 'name', 'B0 map [Hz]')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'B0Map'),'-dtiff')
end

%% %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%              %%%
%%%  STATISTICS  %%%
%%%              %%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATISTICS IN WHOLE SLICES  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Statistics.noBins       =200;
Statistics.B0Data       =B0.Map_Hz(logical(B0.Map_Hz));
Statistics.B0Min        =min(Statistics.B0Data);
Statistics.B0Max        =max(Statistics.B0Data);
Statistics.BinSize_Hz   =(Statistics.B0Max-Statistics.B0Min)/(Statistics.noBins-1);
Statistics.Bins         =Statistics.B0Min :Statistics.BinSize_Hz : Statistics.B0Max;% 1x300 vector
Statistics.B0hist       =histc(Statistics.B0Data,Statistics.Bins);
Statistics.B0hist_norm  =Statistics.B0hist/length(Statistics.B0Data);
Statistics.B0abs_avg    =mean(Statistics.B0Data);
Statistics.B0abs_std    =std(Statistics.B0Data);
Statistics.B0_CoVa        =Statistics.B0abs_std/Statistics.B0abs_avg;


figure(6)
bar(Statistics.Bins,Statistics.B0hist_norm);
axis square
set(gcf, 'name', 'B0+ histogram')
set(gca,'XLim',[Statistics.B0Min Statistics.B0Max])
xlabel('B_0 distribution [Hz]')
ylabel('Normalized Distribution')
set(gca,'Title',text('String',['AVG = ', num2str(Statistics.B0abs_avg), ' Hz , STD = ',num2str(Statistics.B0abs_std),' Hz , bin size = ',num2str(Statistics.BinSize_Hz),' Hz']))
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'B0hist'),'-dtiff')
end

%% %%%%%%%%%%%%%%%%%%%%%
% STATISTICS INSIDE ROI%
%%%%%%%%%%%%%%%%%%%%%%%%
% take only data inside ROI and without outliers
if (Info.enROI)
    
    temp=B0.Map_Hz.*Masks.ROI .* ~Masks.B0outliers;
    temp2=temp(logical(temp));
    Statistics.ROI.B0Data       =temp2;
    Statistics.ROI.B0Max        =max(temp2);
    Statistics.ROI.B0Min        =min(temp2);
    Statistics.ROI.BinSize_Hz   =(Statistics.ROI.B0Max-Statistics.ROI.B0Min)/(Statistics.noBins-1);
    Statistics.ROI.Bins         =Statistics.ROI.B0Min :Statistics.ROI.BinSize_Hz : Statistics.ROI.B0Max;
    Statistics.ROI.B0hist       =histc(temp2,Statistics.ROI.Bins);
    Statistics.ROI.B0hist_norm  =Statistics.ROI.B0hist/length(temp2);
    Statistics.ROI.B0abs_avg    =mean(temp2);
    Statistics.ROI.B0abs_std    =std(temp2);
    Statistics.ROI.B0_CV        =Statistics.ROI.B0abs_std/Statistics.ROI.B0abs_avg;
    clear temp*
    
    figure(7)
    bar(Statistics.ROI.Bins,100*Statistics.ROI.B0hist_norm);
    axis square
    set(gcf, 'name', 'B_0 histogram inside ROI')
    set(gca,'XLim',[Statistics.ROI.B0Min Statistics.ROI.B0Max])
    xlabel('B_0 distribution inside ROI [Hz]')
    ylabel('Normalized Distribution [%]')
    %set(gca,'Title',text('String',['AVG = ', num2str(Statistics.ROI.B0abs_avg),' Hz , STD = ',num2str(Statistics.ROI.B0abs_std),' Hz']))
    set(gca,'Title',text('String',['AVG = ', num2str(Statistics.ROI.B0abs_avg),' Hz , STD = ',num2str(Statistics.ROI.B0abs_std),' Hz , bin size = ',num2str(Statistics.ROI.BinSize_Hz),' Hz']))
    if Info.SavingString
        Image_save(fullfile(Info.ProcessedPath,'B0hist_ROI'),'-dtiff')
    end
end

%% cleanup workspace, prepare export and save data
%clean up workspace
clearvars -except B0 Info Magnitude1 Masks Statistics MrProt
if Info.SavingString
    save(fullfile(Info.ProcessedPath,[Info.SavingString,'_',Info.MID]),'-v7.3')
end

B0Map_Hz=B0.Map_Hz;
detailedB0postprocessingData=struct('B0',B0,'Statistics',Statistics,'MrProt',MrProt,'Magnitude_Echo1',Magnitude1,'Masks',Masks,'Info',Info);

if Info.SavingString
    save(fullfile(Info.ProcessedPath,['B0Maps_',Info.MID]),'B0Map_Hz','-v7.3')
    disp('########################################')
    disp('data and images saved to: ')
    disp(Info.ProcessedPath)
    disp('########################################')
end
disp('===>> POSTPROCESSING FINISHED!')

end
