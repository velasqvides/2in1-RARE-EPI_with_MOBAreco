function [Maps,detailedB1postprocessingData]=B1Mapping_20150521(combinedImage,MrProt,newDataDims,config)
% This function calculates B0, normalized B1+ maps, and absolute B1+ maps for
% various phase based methods such as Bloch-Siegert-Shift, adiabatic Bloch-
% Siegert Shift, PhiFA-Cup and some more
%
% INPUT: 
% combinedImage = complex images with four echoes
% newDataDims = cell array with description of the data dimensions
% MrProt = Struct with MR protocol parameters
% config = struct with information such as filename, pathname and savingString 
%          if no saving string is defined, no data will be saved
%          if config.postprocessing.B1_MAPPING_savingSTRING = 'someString' all postprocessing data will be saved in someString_MIDxx.mat
%          format of config is determined by processSiemensRaw script!
%
% OUTPUT: 
% Maps.B1Map_rel = relative B1 map (wherever B1Map_rel=100% => actuaal FA = nominal FA of protocol)
% Maps.B1Map_uT = absolute B1+ map in [uT/sqrt(kW)]
% Maps.B0Map_Hz = B0 map in [Hz]
% detailedB1postprocessingData = structure of other relevant data
%
% changelog: see end of file
%
% 2014 Matthias Dieringer
% matthias.dieringer@charite.de
%
% modifications by Olli Weinberger
% 2015.05.21
if(nargin~=4)
    disp('Wrong number of arguments!')
    help B1Mapping
    return
end
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

if isempty(config.postprocessing.B1_MAPPING_savingSTRING)
    Info.SavingString=0;
else
    Info.SavingString=genvarname(config.postprocessing.B1_MAPPING_savingSTRING);
    mkdir(Info.ProcessedPath);
end
%% for the B1+ Mapping routine a certain order of the indices is more convenient
CellOrder={'Col','Lin','Eco','Par','Sli'};

% delete 'Channels' from the dimensions because the channels were already combined
newDataDims(strcmp(newDataDims,'Cha'))=[];
%get the permutation indices to get to the new order
for n=1:length(CellOrder)
    CellPos=find(~cellfun(@isempty,strfind(newDataDims,CellOrder{n})));
    if ~isempty(CellPos)
        permOrder(n)=CellPos;
    end
end
permOrder=permOrder(permOrder~=0);
%permute the data according to the new order
CmplxImage=permute(combinedImage,permOrder);
MagMax=max(abs(CmplxImage(:)));
CmplxImage=CmplxImage/MagMax; % normalize magnitude to 1
% extract dimensions into variables
nPar=size(CmplxImage,4); % number of 3D partitions/slices

%% get meta data
[Pulse,Sequence]=getRAWmetadata(MrProt);

% only if B1-Mapping Sequence Version 2013-07-02 was used
%[Pulse,Sequence]=getRAWmetadata_OK(MrProt);
%Pulse.TotalDuration=Pulse.LoopDUR;

%% calculate phase differences for B0 map
% "phasedifference_B0 = phase_Echo2 - phase_Echo1 = phase_Echo4 - phase_Echo3"
% "phasedifference_B1 = phase_Echo1 - phase_Echo3 = phase_Echo2 - phase_Echo4"

%imPhaDiffB0=squeeze(mod(angle(CmplxImage(:,:,2,:))-angle(CmplxImage(:,:,1,:))+pi,2*pi)-pi);
imPhaDiffB0=squeeze(phaseDiff(CmplxImage(:,:,2,:),CmplxImage(:,:,1,:)));
% imPhaDiffB0=squeeze(phaseDiff(CmplxImage(:,:,4,:),CmplxImage(:,:,3,:)));

%imPhaDiffB1=squeeze(mod(angle(CmplxImage(:,:,1,:))-angle(CmplxImage(:,:,3,:))+pi,2*pi)-pi);
imPhaDiffB1=squeeze(phaseDiff(CmplxImage(:,:,1,:),CmplxImage(:,:,3,:)));
%imPhaDiffB1=squeeze(phaseDiff(CmplxImage(:,:,2,:),CmplxImage(:,:,4,:)));

Magnitude1=squeeze(abs(CmplxImage(:,:,1,:)));

%% %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%           %%%
%%%  MASKING  %%%
%%%           %%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% this method is robust even for B1voids but can't cut away inner parts (eg lung tissue)
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
        Masks.lineOfROI=getPosition(ROI);
        %ROI=imellipse(figID);
        %ROI=imrect(figID);
        wait(ROI);
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

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    %%%
%%%  PHASE UNWRAPPING  %%%
%%%                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nPar~=1 % multislice or 3D B1Map
    for cPar=1:nPar % display B0 phase differences to check for phase wraps
        data2plot=imPhaDiffB0(:,:,cPar).*Masks.MagThres(:,:,cPar);
        Min=min(data2plot(:));
        Max=max(data2plot(:));
        Image_plot3D(21,nPar,cPar,Min,Max,'jet',data2plot,0,1)
    end
    for cPar=1:nPar % display B1 phase differences to check for phase wraps
        data2plot=imPhaDiffB1(:,:,cPar).*Masks.MagThres(:,:,cPar);
        Min=min(data2plot(:));
        Max=max(data2plot(:));
        Image_plot3D(22,nPar,cPar,Min,Max,'jet',data2plot,0,1)
    end
    % Show magnitude and phase (does everything look good?) (optional)
    % does your phase data need unwrapping?
    disp('Please check B0 and B1 phase differences for phase wraps!');
    Info.enUnwrapping=input('Does your data need phase unwrapping [1] or not [Enter]? ');
else
    Info.enUnwrapping=1;
end
close all
if (~isempty(Info.enUnwrapping))
    for cPar=1:nPar
        subplot(1,2,1)
        imagescwithnan(imPhaDiffB0(:,:,cPar)./Masks.MagThres(:,:,cPar),jet(256),[1 1 1])
        set(gca,'Title',text('String',['Phase difference used for B_0-Map. slice: ', num2str(cPar)]));
        axis image
        subplot(1,2,2)
        imagescwithnan(imPhaDiffB1(:,:,cPar)./Masks.MagThres(:,:,cPar),jet(256),[1 1 1])
        set(gca,'Title',text('String','Phase difference used for B_1-Map'))
        axis image
        colormap(jet)
        inpt=input('Are phase wraps in B0 [1], B1 [2], both [3] or none [Enter]? ');
        
        % perform unwrapping
        if (~isempty(inpt))
            cmplx1=Magnitude1(:,:,cPar).*exp(1i.*imPhaDiffB0(:,:,cPar));
            cmplx2=Magnitude1(:,:,cPar).*exp(1i.*imPhaDiffB1(:,:,cPar));
            
            while 1
                if (inpt==1 || inpt==3)
                    close all
                    imagesc(angle(cmplx1))
                    colorbar
                    axis image
                    set(gca,'Title',text('String','Pick a point, where \Delta B_0 is small and with good signal!'))
                    [xpoint,ypoint] = ginput(1);
                    disp('===>> Unwrapping, please wait...')
                    imPhaDiffB0(:,:,cPar)=Unwrap2D_QualityGuided_inputXY(cmplx1, Info.relativeMagnitudeThreshold ,0,xpoint,ypoint);
                    disp('===>> Unwrapping done!')
                end
                
                if (inpt==2 || inpt==3)
                    close all
                    imagesc(angle(cmplx2))
                    colorbar
                    axis image
                    set(gca,'Title',text('String','Pick a point, where B_1^+ is small but sufficient signal!'))
                    [xpoint,ypoint] = ginput(1);
                    disp('===>> Unwrapping, please wait...')
                    imPhaDiffB1(:,:,cPar)=Unwrap2D_QualityGuided_inputXY(cmplx2, Info.relativeMagnitudeThreshold ,0,xpoint,ypoint);
                end
                
                subplot(1,2,1)
                imagesc(imPhaDiffB0(:,:,cPar))
                colormap jet
                set(gca,'Title',text('String',['Phase difference used for B_0-Map, slice: ' num2str(cPar)]))
                axis image
                colorbar
                subplot(1,2,2)
                imagesc(imPhaDiffB1(:,:,cPar))
                set(gca,'Title',text('String','Phase difference used for B_1-Map'))
                axis image
                colorbar
                
                inpt2=input('DONE! Happy with the result [any key] or not [1]? ');
                %inpt2=[];
                
                inpt2(isempty(inpt2))=0;
                if inpt2==0
                    break
                end
            end
        end
    end
end
Info.UnwrappingInput=inpt;
%% prepare for export [deg]
close all
% cleanup workspace
clearvars -except Info Masks Magnitude1 Pulse Sequence nPar imPhaDiffB0 imPhaDiffB1 MrProt

%Display meta data
Pulse
Sequence

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  B0 CALCULATIONS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate B0 Map
disp('===>> Calculating B0 maps...');
B0.Map=imPhaDiffB0/2/pi/Sequence.DeltaTE .*Masks.MagThres;
Info.MaxDeltaB0=max(abs(B0.Map(:)));
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-Info.MaxDeltaB0,Info.MaxDeltaB0,jet,B0.Map(:,:,cPar)./Masks.MagThres(:,:,cPar),1,1)
end
set(gcf, 'name', 'B0 map')

%cut off outliers
disp(['Maximum B0 is ',num2str(round(Info.MaxDeltaB0)),' Hz']);
Info.B0cut=input('choose cut off value ([Enter] to use max)> ');
Info.B0cut(isempty(Info.B0cut))=Info.MaxDeltaB0;
Info.B0cut=Info.MaxDeltaB0;
B0.outliers=logical(abs(B0.Map)>Info.B0cut);
B0.noOutliers=sum(B0.outliers(:));
B0.Map(abs(B0.Map)>Info.B0cut)=0;
Masks.B0outliers=B0.outliers;

%% save image and data
close gcf
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-Info.B0cut,Info.B0cut,'jet',B0.Map(:,:,cPar).*(1-B0.outliers(:,:,cPar)),1,1)
end
set(gcf, 'name', 'B0 map [Hz]')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'B0Map'),'-dtiff')
end
close (gcf)
%% %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  B1 CALCULATIONS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===>> Calculating B1+ maps...');
%B1.Sub=imPhaDiffB1*180/pi .*Masks.MagThres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  1)Sensitivity curve  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%round up B0 distribution and add +20Hz
SimulatedSensitivity.DistrB0=ceil(max(abs(B0.Map(:)))/10)*10;

%number of spins adapted to data
SimulatedSensitivity.SpinsB0=SimulatedSensitivity.DistrB0/10;% this is good for Bloch Siegert
SimulatedSensitivity.B1Max=4; % 400% is a realistic upper limit. but more can be feasible for some surface coils!
SimulatedSensitivity.SpinsB1=100*SimulatedSensitivity.B1Max;%number of spins adapted
Pulse.Shift=0;

%optional: relaxation times of tissue. They only have a minor effect on bloch simulations. 
SimulatedSensitivity.T1=30e-3;%[s]
SimulatedSensitivity.T2=1e-3; %[s]

% calculate sensitivity curve with mex function -> takes some time!
SimulatedSensitivity=Sensitivity_build_mex(Pulse,SimulatedSensitivity);

% sometimes due to wrong unwrapping the sensitivity curve starts at phase accrual of 360° and not at 0°.Then you have to manually subtract 360°!
%SimulatedSensitivity.DeltaPhi=SimulatedSensitivity.DeltaPhi-360;


SimulatedSensitivity.B0Range=[-SimulatedSensitivity.DistrB0 +SimulatedSensitivity.DistrB0];
SimulatedSensitivity.B1Range=[0 SimulatedSensitivity.B1Max];
SimulatedSensitivity.DeltaPhiRange=[min(SimulatedSensitivity.DeltaPhi(:)) max(SimulatedSensitivity.DeltaPhi(:))];

%%{
% plot and save sensitivity curve
% if Info.SavingString
%     save(fullfile(Info.ProcessedPath,'Sensitivity'),'SimulatedSensitivity','-v7.3')
% end
figure(3)
surf(SimulatedSensitivity.B0mat,SimulatedSensitivity.B1mat,SimulatedSensitivity.DeltaPhi)
axis square
set(gca,'xLim',[-SimulatedSensitivity.DistrB0-50 SimulatedSensitivity.DistrB0+50],'yLim',[0 max(SimulatedSensitivity.B1mat(:))],'zLim',[min(SimulatedSensitivity.DeltaPhi(:))-10 max(SimulatedSensitivity.DeltaPhi(:))+10],'View',[-35 25],'FontSize',16)
zlabel('Phase accrual \psi [deg]'); ylabel('B_1 [normalized]'); xlabel('B_0 offset [Hz]')
set(gcf, 'name', 'B1+ sensitivity map: Phase accrual')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'Sensitivity'),'-dtiff')
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  2) Interpolation  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for cPar=1:nPar
    %B1map=Interpolation_B1B0(B0map,B1dif,Curve,B0vec,B1vec)
    %B1map                                      =  Interpolation_B1B0(B0map,           B1dif,                                                  Curve,                           B0vec,                          B1vec)
    B1.Map_rel(:,:,cPar)=Masks.MagThres(:,:,cPar).*Interpolation_B1B0(B0.Map(:,:,cPar),imPhaDiffB1(:,:,cPar)*180/pi .*Masks.MagThres(:,:,cPar),...
                                                                       SimulatedSensitivity.DeltaPhi,SimulatedSensitivity.B0vec,SimulatedSensitivity.B1vec);
end
B1.Map_rel(B0.outliers)=NaN;
B1.Singularities= ~isfinite(B1.Map_rel);
B1.Map_rel=setNAN_Inf_to_0(B1.Map_rel);%get rid of NaNs
B1.noSingularities=sum(B1.Singularities(:));
Masks.B1Sing=B1.Singularities;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate absolute B1+ map in uT/sqrt(kW) & Masks.B0_B1_Mag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B1.rel2absFactor = Pulse.NominalB1 / Sequence.ReqVolt * sqrt(50) * sqrt(1000);
% = const * (number of half-loops)*(duration of one loop)/(duration of entire pulse)/(reference amplitude)
% independent of loopFA !!!
B1.Map_uT = B1.Map_rel * B1.rel2absFactor;%unit: uT/sqrt(kW)
for cPar=1:nPar
    Image_plot3D(4,nPar,cPar,0,3,'hot',B1.Map_rel(:,:,cPar),1,0)
end
set(gcf, 'name', 'normalized B1+ [%/100]')
%%
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'B1Map_rel'),'-dtiff')
end
Info.MaxB1_uT=max(B1.Map_uT(:));
for cPar=1:nPar
    Image_plot3D(5,nPar,cPar,0,Info.MaxB1_uT,'hot',B1.Map_uT(:,:,cPar),1,0)
    if (Info.enROI)
        hold on
        ROIcontour=plot(Masks.lineOfROI(:,1),Masks.lineOfROI(:,2));
        set(ROIcontour,'Color','blue','LineWidth',2);
    end
end
set(gcf, 'name', 'Absolute B1+ [uT/sqrt(kW)]')

if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'B1Map_uT'),'-dtiff')
end

%% %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%              %%%
%%%  STATISTICS  %%%
%%%              %%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

% combining all Masks
Masks.B0_B1_Mag= ~Masks.B1Sing & Masks.MagThres & ~Masks.B0outliers;
if (Info.enROI)
    Masks.B0_B1_Mag_ROI = Masks.B0_B1_Mag & Masks.ROI;
end
Info.screenPrintout{1,1}=sprintf('Mean normalized B1_wholeSlice = %1.2f +- %1.2f %%',100*mean2(B1.Map_rel(Masks.B0_B1_Mag)),100*std2(B1.Map_rel(Masks.B0_B1_Mag)) );
disp(Info.screenPrintout{1,1})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATISTICS IN WHOLE SLICES  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Statistics.B1Data       =B1.Map_uT(Masks.B0_B1_Mag);
Statistics.B1Max        =max(Statistics.B1Data);
Statistics.Bins         =0 :Statistics.B1Max/(300-1) : Statistics.B1Max;% 1x300 vector
Statistics.B1hist       =histc(Statistics.B1Data,Statistics.Bins);
Statistics.B1hist_norm  =Statistics.B1hist/length(Statistics.B1Data);
Statistics.B1abs_avg    =mean(Statistics.B1Data);
Statistics.B1abs_std    =std(Statistics.B1Data);
Statistics.B1_CV        =Statistics.B1abs_std/Statistics.B1abs_avg;

if ~(Info.enROI)
    figure(6)
    bar(Statistics.Bins(2:end),Statistics.B1hist_norm(2:end));
    axis square
    set(gcf, 'name', 'B1+ histogram')
    set(gca,'XLim',[0 Statistics.B1Max])
    xlabel('B_1^+ distribution [uT]')
    ylabel('Normalized Distribution')
    set(gca,'Title',text('String',['AVG = ', num2str(Statistics.B1abs_avg), ' uT , STD = ' num2str(Statistics.B1abs_std),' uT']))
    if Info.SavingString
        Image_save(fullfile(Info.ProcessedPath,'B1hist'),'-dtiff')
    end
end
%% %%%%%%%%%%%%%%%%%%%%%
% STATISTICS INSIDE ROI%
%%%%%%%%%%%%%%%%%%%%%%%%
% take only data inside ROI & above threshold & no B1 singularity
if (Info.enROI)
    Statistics.ROI.Mask         =Masks.B0_B1_Mag_ROI;
    Statistics.ROI.B1Data       =B1.Map_uT(Masks.B0_B1_Mag_ROI);
    Statistics.ROI.B1Max        =max(Statistics.ROI.B1Data);
    Statistics.ROI.Bins         =(0 : 4/numel(Statistics.ROI.B1Data) :1)* Statistics.ROI.B1Max;% adapt number of bins to available data
    Statistics.ROI.B1hist       =histc(Statistics.ROI.B1Data,Statistics.ROI.Bins);
    Statistics.ROI.B1hist_norm  =Statistics.ROI.B1hist/length(Statistics.ROI.B1Data);
    Statistics.ROI.B1abs_avg    =mean(Statistics.ROI.B1Data);
    Statistics.ROI.B1abs_std    =std(Statistics.ROI.B1Data);
    Statistics.ROI.B1_CV        =Statistics.ROI.B1abs_std/Statistics.ROI.B1abs_avg;
    
    Statistics.ROI.B1rel_avg    =mean(B1.Map_rel(Masks.B0_B1_Mag_ROI));
    
    figure(7)
    bar(Statistics.ROI.Bins(2:end),100*Statistics.ROI.B1hist_norm(2:end));
    axis square
    set(gcf, 'name', 'B1+ histogram inside ROI')
    set(gca,'XLim',[0 Statistics.ROI.B1Max])
    xlabel('B_1^+ distribution inside ROI [uT]')
    ylabel('Normalized Distribution [%]')
    set(gca,'Title',text('String',['AVG = ', num2str(Statistics.ROI.B1abs_avg),' uT, STD = ' num2str(Statistics.ROI.B1abs_std),' uT']))
    if Info.SavingString
        Image_save(fullfile(Info.ProcessedPath,'B1hist_ROI'),'-dtiff')
    end
    Info.screenPrintout{2,1}=sprintf('B1_ROI = %1.2f +- %1.2f uT',Statistics.ROI.B1abs_avg,Statistics.ROI.B1abs_std);
    Info.screenPrintout{3,1}=sprintf('Old Reference Voltage: %3.1f V ',Sequence.RefVolt);
    Info.screenPrintout{4,1}=sprintf('calibrated transmitter voltage would be %3.1f V ',Sequence.RefVolt/Statistics.ROI.B1rel_avg);
    Info.screenPrintout{5,1}=sprintf('Statistics.ROI.B1rel_avg: %1.2f percent ',Statistics.ROI.B1rel_avg);
    disp(Info.screenPrintout{2,1})
    disp(Info.screenPrintout{3,1})
    disp(Info.screenPrintout{4,1})
    disp(Info.screenPrintout{5,1})
end

%% cleanup workspace, prepare export and save data
%clean up workspace
B1.phaseDiff_unwrapped=imPhaDiffB1;
B0.phaseDiff_unwrapped=imPhaDiffB0;
clearvars -except B0 B1 Info Magnitude1 Masks Pulse Sequence SimulatedSensitivity Statistics MrProt

if (Info.enROI)
    Maps.B1Map_rel=B1.Map_rel .*Masks.B0_B1_Mag_ROI;
    Maps.B0Map_Hz=B0.Map .*Masks.B0_B1_Mag_ROI;
    Maps.B1Map_uT=B1.Map_uT .*Masks.B0_B1_Mag_ROI;
else
    Maps.B1Map_rel=B1.Map_rel .*Masks.B0_B1_Mag;
    Maps.B0Map_Hz=B0.Map .*Masks.B0_B1_Mag;
    Maps.B1Map_uT=B1.Map_uT .*Masks.B0_B1_Mag;  
end

if Info.SavingString
    save(fullfile(Info.ProcessedPath,[Info.SavingString,'_',Info.MID]),'-v7.3')
    %save(fullfile('O:\M\projects\ok_4xLocal vs Global_3T\invivoB1-Maps\B1DataPerScan',[Info.SavingString,'_',Info.MID]),'-v7.3')
end

detailedB1postprocessingData=struct('Info',Info,'MrProt',MrProt,'Pulse',Pulse,'Sequence',Sequence,...
    'SimulatedSensitivity',SimulatedSensitivity,'Magnitude_Echo1',Magnitude1,...
    'Masks',Masks,'B0',B0,'B1',B1,'Statistics',Statistics,'Maps',Maps);

if Info.SavingString
    save(fullfile(Info.ProcessedPath,['B0B1Maps_',Info.MID]),'Maps','-v7.3')
    disp('########################################')
    disp('data and images saved to: ')
    disp(Info.ProcessedPath)
    disp('########################################')
end
disp('===>> POSTPROCESSING FINISHED!')
%{
changelog:

2013-05-16 -Adapted from Complete_F_B1new.m to work with Siemens RAW data
            (Matthias Dieringer)
2013-05-22 -added averages
            (Matthias Dieringer)
2013-05-21 -Automatic image registration implemented
            (Matthias Dieringer)
2013-05-29 -changed order of cells to calculate everything with one click
           -Flaw with registration fixed
            (Matthias Dieringer)
2013-05-30 -added phase unwrapping
           -Removed Limit of interpolation curve (AGrässl, Matthias Dieringer)
            that cut off high B1 values.
           -Automated the offset of the interpolation curve. Factors of
            360° are added automatically now (Matthias Dieringer)
           -dramatic speed increase for calculation of the interpolation
            curve. If you don't have the correct MEX-file or can't
            compile it, the calculation will switch back to the slower
            Matlab version and ask for desired B1 scale factor
            (Matthias Dieringer)
2013-06-06 -Speed-up and simplification in the Interpolation_B0B1
function. There were sometimes problems with the lookup
            table (leading to extremely high B1 in some regions),
            which should be fixed now.
            (Matthias Dieringer)
2013-06-13 -changed raw reading routine
           -extended raw data filtering options
           -code is faster and simplified
            (Matthias Dieringer)
2013-06-13 -adapted code to handle 1 receiver only as well
            (Matthias Dieringer)
2013-06-17 -RF pulse and sequence info is fetched from RAW data
            header (Matthias Dieringer)
2013-06-17 -removed some unnecessary code and graphs
           -changed max B1 sensitivity (curve) to 4 with 120 samples
           -changed how the B0 Map is cut off (was with sign switch
            before, now outliers are set to 0)
            (Matthias Dieringer)
2013-07-02 -speed improvements and simplifications in sensitivity_build
            and pulse_build
           -non-MEX version of sensitivity build is not supported anymore
            in this version
           -Adiabatic Bloch Siegert pulse implemented (Khaligi et al, MRM
            2013)
           -Gaussian RF pulse was removed, no application for this
           (Matthias Dieringer)

           THEREFORE YOU NEED TO USE THE NEW VERSION OF THE SEQUENCE
           (20130702 and later)

2013-08-29 -there was a bug in the Pulse_build.m so that the scaling was
            not correct (only the new version 2013-07-02 is affected).
            This could lead to an B1 error of up to 5%
            (Matthias Dieringer)
2013-09-03 -change save(...) to save(..,'-v7.3') to allow huge files and
            to allow loading of single variables out of a saved array
           -change strcat to fullfile in order to make it working for
            Linux as well
           -B0 distribution to simulate is rounded up instead of being
            rounded down (caused NaNs in the interpolation sonmetimes)
            (Olli Kraus)
2013-09-06 -replaced single precision variables by double precision
            variables to make this script compatible for use with Linux
           -Linux 64bit mex bloch simulation file added
2013-09-25 -Added absolute Map in ut/(sqrt/kW). Reference point for the
            RMS input power is the siemens coil plug(A Graessl)
2013-09-27 -meas ID is added to the name of the subfolder "Processed" and
            maps are stored in the "processed" subfolder.
           -after unwrapping, "any key" is used to proceed instead of "1"
           -new function "Display_Image" to display images
           -cleaned up code
           (Matthias Dieringer)
2014-01-13 -3D and multislice is enabled again
           -3D methods using composite RF hard pulses are facilitated now
           -RAW data import and handling adapted
           -MR Sequence was adapted as well, please use the new version
           (although older RAW files should still work)
           -cleaned up some of the code (also external functions)
2014-01-29 -optional ROI drawing with Statistics for that ROI.
           -tested for multi-slice, single slice and 3D (PhiFA-cup)
           -some minor changes
           -all data and images saved in one folder: Processed_MIDxxx_date
           (Olli Kraus)
2014-02-21 -ROI statistics now apply to absolute B1map inside ROI instead
            of relative Map.
           -new parameter in B1-struct: B1.rel2absFactor
           -title('...') replaced by more stable set(gca,'Title',text('String','...'));
            (Olli Kraus)
2014-02-28 -statistics: adapt number of bins to number of available data points
           -speed up by avoiding unnecessary input questions
            (Olli Kraus)
2014-03-04 -struct Masks changed
           -new substrucure added: Statistics.ROI.B1rel_avg -> helpful for
            calculation of FAact of other sequences
           -more information in Info-struct
            (Olli Kraus)
2014-04-15 -screen printouts are now also saved in Info.screenPrintout-struct
            (Olli Kraus)
2014-05-06 -no Sensitivity.tif will be saved
            -enSaving (logical) is now replaced by SavingString (string). Empty string means no saving
            - allData.mat now is renamed to SavingString_MIDxxx.mat
            - tested with 2D raw data. There might be problems with 3D or multislice data
            (Olli Kraus)
2014-09-16  -optional input argument config.postprocessing.B1_MAPPING_savingSTRING:
                if ='someString' -> data and images are beeing saved
                if =[]           -> no saving
2016-09-21 -now subfunctions also support sodium imaging:
            changes in bloch.c/bloch.mexw64 and Pulse_build.m were necessary
            (Olli Weinberger)
%}
end
