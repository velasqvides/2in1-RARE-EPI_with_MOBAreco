function Maps=B1Mapping_DAM_45_90(Seq1,Seq2,SavingString)
% changelog: see end of file
%
% written by Olli Weinberger
% 2016.02.11

%%
Info.usedMfile=[mfilename('fullpath'),'.m'];
Info.TimeOfExecution=datestr(clock);

Info.usedSeq=Seq1.MrProt.SequenceFileName;

Info.RawPath=Seq1.config.pathname;
temp1=regexpi(Seq1.config.filename,'_','split');
Info.MID1=temp1{2};

temp2=regexpi(Seq2.config.filename,'_','split');
Info.MID2=temp2{2};

Info.MIDs=[Info.MID1,'_',Info.MID2];

Info.SavingString=SavingString;
Info.ProcessedFolder=strcat(Info.MIDs ,'_Processed_',datestr(now,'yyyymmdd'));
Info.ProcessedPath=fullfile(Seq1.config.pathname,Info.ProcessedFolder);
if Info.SavingString
    mkdir(Info.ProcessedPath);
end
% if field is not defined: set default value (->'calculate')
if not(isfield(Seq1.config.postprocessing,'B1_Mapping_RFpulseShape'))
    Seq1.config.postprocessing.B1_Mapping_RFpulseShape = 'calculate';
end
if not(isfield(Seq2.config.postprocessing,'B1_Mapping_RFpulseShape'))
    Seq2.config.postprocessing.B1_Mapping_RFpulseShape = 'calculate';
end
%% determine which sequence is with FA and which with 2*FA
% 2*FA_1=FA_2
temp1=Seq1.MrProt.FlipAngleDegree;
temp2=Seq2.MrProt.FlipAngleDegree;
if temp2==2*temp1
    Seq_1=Seq1;
    Seq_2=Seq2;
elseif 2*temp2==temp1
    Seq_2=Seq1;
    Seq_1=Seq2;
else
    disp('wrong FAs of protocols')
end
%Info.reconData_Seq1 = Seq_1.Seq1reconData;
%Info.reconData_Seq2 = Seq_2.Seq1reconData;

%=> now, Seq_1 is FA and Seq_2 is with 2*FA
clear Seq1 Seq2 temp*
%% read out important sequence parameters
% for sequence with FA1

FA_1=Seq_1.MrProt.FlipAngleDegree;
RefAmp_1=Seq_1.MrProt.TXSPEC.NucleusInfo(1,1).ReferenceAmplitude;
RealAmp_1=Seq_1.MrProt.TXSPEC.RFPULSE.Amplitude;
Mag_1=cast(abs(Seq_1.combinedImage),'double');
Pha_1=cast(angle(Seq_1.combinedImage),'double');

FA_2=Seq_2.MrProt.FlipAngleDegree;
RefAmp_2=Seq_2.MrProt.TXSPEC.NucleusInfo(1,1).ReferenceAmplitude;
RealAmp_2=Seq_2.MrProt.TXSPEC.RFPULSE.Amplitude;
Mag_2=cast(abs(Seq_2.combinedImage),'double');
Pha_2=cast(angle(Seq_2.combinedImage),'double');

%% check if parameters are consistent

if RefAmp_1~=RefAmp_2
    disp('Reference Amplitudes differ!')
end

if 2*RealAmp_1~=RealAmp_2
    disp('Pulse-Amplitudes were not scaled linearly!')
end
disp(['RealAmp_1 = ',int2str(RealAmp_1),' V'])
disp(['RealAmp_2 = ',int2str(RealAmp_2),' V'])

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
[Masks.MagThres,Masks.Image,Info.relativeMagnitudeThreshold]=masking_by_thresholding_3D(Mag_1);
Masks.MagThres=logical(Masks.MagThres);
close(gcf)
nPar=Seq_1.twix_obj.image.NPar;                                             %If twix emty --> nPar = 1
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
%Info.enROI=[];
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
        %wait(ROI);
        Masks.ROI(:,:,cPar)=createMask(ROI);
        Image_plot3D(2,nPar,cPar,0,1,'jet',Masks.ROI(:,:,cPar).* Masks.Image(:,:,cPar),0,1)
    end
    set(gcf, 'name', 'ROIMask')
    if Info.SavingString
        Image_save(fullfile(Info.ProcessedPath,'ROIMask'),'-dtiff')
    end
else
    Info.enROI=0;
    disp('no ROI drawn')
    Masks.ROI=true(size(Masks.MagThres));
end
close all

%% read out the relevant parameter from MrProt
nuc_1=Seq_1.MrProt.TXSPEC.NucleusInfo(1,1).Nucleus;
nuc_2=Seq_2.MrProt.TXSPEC.NucleusInfo(1,1).Nucleus;

GradientMode_1=Seq_1.MrProt.GRADSPEC.Mode;
GradientMode_2=Seq_2.MrProt.GRADSPEC.Mode;

RFMode_1=Seq_1.MrProt.TXSPEC.RFPulseType;
RFMode_2=Seq_2.MrProt.TXSPEC.RFPulseType;

if isequal(GradientMode_1,GradientMode_2,1)
    GradientMode='fast';
elseif isequal(GradientMode_1,GradientMode_2,2)
    GradientMode='normal';
elseif isequal(GradientMode_1,GradientMode_2,4)
    GradientMode='whisper';
else
    disp('Different gradient modes were used for Seq1 and Seq2!')
    return
end

if isequal(nuc_1,nuc_2,'23Na')
    gamma=7.09051647E7; %[rad/s/T]
elseif isequal(nuc_1,nuc_2,'1H')
    gamma=2.67522212E8; %[rad/s/T]
else
    disp('some error with the used nuclei happened!')
    return
end
clear GradientMode_1 GradientMode_2 
disp(['used nuclei: ',nuc_1])

%% calculate the properties of the RF pulse
disp([Seq_1.config.postprocessing.B1_Mapping_RFpulseShape, ' RF pulse'])

if strcmp(Seq_1.config.postprocessing.B1_Mapping_RFpulseShape,'measured') || strcmp(Seq_2.config.postprocessing.B1_Mapping_RFpulseShape,'measured')
    RFpulseCharacteristics=calculatePulseProperties_fromMeasuredPulse(Seq_1.config.postprocessing.B1_Mapping_RFpulse_filename,RealAmp_1);

elseif strcmp(Seq_1.config.postprocessing.B1_Mapping_RFpulseShape,'simulated') || strcmp(Seq_2.config.postprocessing.B1_Mapping_RFpulseShape,'simulated')

    RFpulseCharacteristics=calculatePulseProperties_fromSimulatedPulse(Seq_1.config.postprocessing.B1_Mapping_RFpulse_filename,RealAmp_1);
    
elseif strcmp(Seq_1.config.postprocessing.B1_Mapping_RFpulseShape,'calculate') || strcmp(Seq_2.config.postprocessing.B1_Mapping_RFpulseShape,'calculate')
    
    if isequal(RFMode_1,RFMode_2,1)
        temp_RFMode='fast'
        %temp_tau_pulse_s=400e-6; % if BWrx>1560Hz/pixel
        temp_tau_pulse_s=1000e-6;% if BWrx<1560Hz/pixel
        temp_BWT=2; %BWT=transmit bandwith-time-product
    elseif isequal(RFMode_1,RFMode_2,2)                                     %Pay attention to the Pulse Duration
        temp_RFMode='normal'
        temp_tau_pulse_s=600e-6;
%         temp_tau_pulse_s=2000e-6;
        temp_BWT=2.7;
    elseif isequal(RFMode_1,RFMode_2,4)
        temp_RFMode='lowSAR'
        temp_tau_pulse_s=4e-3;
        temp_BWT=2;
    else
        disp('Different RF modes were used for Seq1 and Seq2!')
        return
    end
    RFpulseCharacteristics=calculatePulseProperties_HanSinc_BWT(RealAmp_1,temp_tau_pulse_s,temp_BWT);
    RFpulseCharacteristics.RFMode=temp_RFMode;
else
    disp('unknown RF pulse shape')
    disp(Seq_1.MrProt.TXSPEC.RFPULSE.Name)
end
clear temp* RFMode_1 RFMode_2
%% %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  FLIP ANGLE MAPS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===>> Calculating B1+ maps...');
%B1.Sub=imPhaDiffB1*180/pi .*Masks.MagThres;

%%% THEORY:
%%%
%%% Mag1(FA1)/Mag2(FA2=2*FA1) = 1/2/cos(FAmap)
%%%
%%%  0?<=FA<=90? => (1)
%%% ==> FAmap= acos (Mag2/2/Mag1)
%%%
%%%  90?<=FA<=180? => (2)
%%% ==> FAmap= pi - acos (Mag2/2/Mag1)

%% define regions where (1) or (2) or none of these
tolAngle=45;

Diff_Pha=180/pi*phaseDiff(Seq_2.combinedImage,Seq_1.combinedImage);%uses atan2-function
%  0?<=FA<=90? => (1)
Masks.tolAngle=tolAngle;
Masks.phaseRange_0_90=false(size(Masks.Image));
Masks.phaseRange_0_90(Diff_Pha <= tolAngle & -tolAngle<=Diff_Pha)=true;
%  90?<=FA<=180? => (2)
Masks.phaseRange_90_180=false(size(Masks.Image));
Masks.phaseRange_90_180(Diff_Pha>=180-tolAngle | Diff_Pha<=tolAngle-180)=true;
% out of tolerance region
Masks.phaseRange_offPhaseTol=zeros(size(Masks.Image));
Masks.phaseRange_offPhaseTol=~(Masks.phaseRange_0_90 | Masks.phaseRange_90_180);


%% actual FA calculation
temp1=acos((Masks.phaseRange_0_90.*Mag_2) ./(2*(Masks.phaseRange_0_90.*Mag_1)));        %   0?<=FA<=90?  => (1)
temp2=pi-acos((Masks.phaseRange_90_180.*Mag_2) ./(2*(Masks.phaseRange_90_180.*Mag_1))); %  90?<=FA<=180? => (2)
temp1=setNAN_Inf_to_0(temp1);
temp2=setNAN_Inf_to_0(temp2);
B1.FA_Map_rad=real(temp1)+real(temp2);%temp1 and temp2 are complementary
B1.FA_Map_degree=180/pi.*B1.FA_Map_rad;
B1.FA_Map_rel = B1.FA_Map_degree ./ FA_1; %relative FA-map.
%where FA_rel=1 nominal flip angle of protocol is actually achieved flip angle.
B1.Singularities=logical(imag(temp1) | imag(temp2));
B1.noSingularities=sum(B1.Singularities(:));
temp=B1.Singularities .* Masks.ROI;
B1.noSingularitiesROI=sum(temp(:));
Masks.B1Singularities=~B1.Singularities;

%% combining all Masks
%Masks.B1_Mag_Phase = ~Masks.phaseRange_offPhaseTol & Masks.B1Singularities & Masks.MagThres;
Masks.B1_Mag_Phase = Masks.phaseRange_0_90 & Masks.B1Singularities & Masks.MagThres;
Masks.B1_Mag_Phase_ROI = Masks.B1_Mag_Phase & Masks.ROI;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               %%%
%%%  CALCULATE ABSOLUTE B1+ MAPS  %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate PSI: conversion factor from currents in the RF coil to produced B1+ field in the sample
% PSI=PSI(r) carries all the spatial dependence (but not the temporal dependence)
% unit of PSI: nT/V
%
% FA_Map_rad (r) = gamma    * integral_0..tau (B1_uT(r,t)dt) =
%                = gamma    * psi (r) * integral_0..tau (Waveform_V(t)dt)
% Units:
%  [rad]         = [rad/s/T] * [T/V] * [V*s]= [rad]

B1.offPhaseTol=Masks.phaseRange_offPhaseTol;
B1.no_offPhaseTol=sum(Masks.phaseRange_offPhaseTol(:));
B1.B1_nT_V  = B1.FA_Map_rad /gamma / RFpulseCharacteristics.pulseIntegral_Vs *1e9;
% 1kW @ 50 Ohm <=> sqrt(1000W*50Ohm)=2*sqrt(5)*V=224V
% x uT/sqrt(kW) = x*2*sqrt(5) nT/V
B1.B1_uT_kW = B1.B1_nT_V /2/ sqrt(5);

%% plot maps and save as tiffs
Info.MaxB1_uT_kW=max(B1.B1_uT_kW(:));
Info.MaxFA_degree=max(B1.FA_Map_degree(:));
% relative FAmap
for cPar=1:nPar
    Image_plot3D(3,nPar,cPar,0,3,'hot',B1.FA_Map_rel(:,:,cPar).*Masks.B1_Mag_Phase_ROI,1,0)
end
set(gcf, 'name', 'relative FA-Map [100%]')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'relativeFA_Map'),'-dtiff')
end
% absolute FAmap
for cPar=1:nPar
    Image_plot3D(4,nPar,cPar,0,Info.MaxFA_degree,'hot',B1.FA_Map_degree(:,:,cPar).*Masks.B1_Mag_Phase_ROI,1,0)
end
set(gcf, 'name', 'FA-Map [?]')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'FA_Map_degree'),'-dtiff')
end
% absolute B1map
for cPar=1:nPar
    Image_plot3D(5,nPar,cPar,0,Info.MaxB1_uT_kW,'hot',B1.B1_uT_kW(:,:,cPar).*Masks.B1_Mag_Phase_ROI,1,0)
    if (Info.enROI)
        hold on
        ROIcontour=plot(Masks.lineOfROI(:,1),Masks.lineOfROI(:,2));
        set(ROIcontour,'Color','black','LineWidth',2);
    end
end
set(gcf, 'name', 'Absolute B1+ [uT/sqrt(kW)]')
if Info.SavingString
    Image_save(fullfile(Info.ProcessedPath,'B1Map_uT_kW'),'-dtiff')
    Image_save(fullfile(Info.ProcessedPath,'B1Map_uT_kW'),'fig')
    %Image_save(['O:\M\projects\ok_4xLocal vs Global_3T\invivoB1-Maps\Tifs\',Info.SavingString,'_',Info.MID],'-dtiff')
end
Info.screenPrintout{1,1}=sprintf('Mean B1+ = %1.2f +- %1.2f uT/sqrt(kW) ',mean2(B1.B1_uT_kW(Masks.B1_Mag_Phase)),std2(B1.B1_uT_kW(Masks.B1_Mag_Phase)) );
disp(Info.screenPrintout{1,1})

%% %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%              %%%
%%%  STATISTICS  %%%
%%%              %%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
Statistics.B1Data       =B1.B1_uT_kW(Masks.B1_Mag_Phase);
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
    Statistics.ROI.B1Data       =B1.B1_uT_kW(Masks.B1_Mag_Phase_ROI);
    Statistics.ROI.B1Max        =max(Statistics.ROI.B1Data);
    Statistics.ROI.Bins         =(0 : 4/numel(Statistics.ROI.B1Data) :1)* Statistics.ROI.B1Max;% adapt number of bins to available data
    Statistics.ROI.B1hist       =histc(Statistics.ROI.B1Data,Statistics.ROI.Bins);
    Statistics.ROI.B1hist_norm  =Statistics.ROI.B1hist/length(Statistics.ROI.B1Data);
    Statistics.ROI.B1abs_avg    =mean(Statistics.ROI.B1Data);
    Statistics.ROI.B1abs_std    =std(Statistics.ROI.B1Data);
    Statistics.ROI.B1_CV        =Statistics.ROI.B1abs_std/Statistics.ROI.B1abs_avg;
    
    Statistics.ROI.B1rel_avg    =mean(B1.FA_Map_rel(Masks.B1_Mag_Phase_ROI));
    
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
    Info.screenPrintout{2,1}=sprintf('(B1act/B1nom)_ROI = %1.2f +- %1.2f uT',Statistics.ROI.B1abs_avg,Statistics.ROI.B1abs_std);
    Info.screenPrintout{3,1}=sprintf('Old Reference Voltage: %3.1f V ',RefAmp_1);
    Info.screenPrintout{4,1}=sprintf('Adjust Reference Voltage to %3.1f V ',RefAmp_1/Statistics.ROI.B1rel_avg);
    disp(Info.screenPrintout{2,1})
    disp(Info.screenPrintout{3,1})
    disp(Info.screenPrintout{4,1})
end
%% cleanup workspace, prepare export and save data
Seq1=Seq_1;
Seq2=Seq_2;
clear temp* ROI* cPar nPar BWT *Mode *_1 *_2

% structure for direct output
Maps.FA_rel=B1.FA_Map_rel.*Masks.B1_Mag_Phase_ROI;
Maps.FA_degree=B1.FA_Map_degree.*Masks.B1_Mag_Phase_ROI;
Maps.B1_uT_kW=B1.B1_uT_kW.*Masks.B1_Mag_Phase_ROI;

% if Info.SavingString
%     save(fullfile(Info.ProcessedPath,[Info.SavingString,'_',Info.MIDs]),'-v7.3')
% end
% 
% if Info.SavingString
%     save(fullfile(Info.ProcessedPath,['FA_B1_Maps_',Info.MIDs]),'Maps','-v7.3')
% end
% 
% disp('data and images saved to: ')
% disp(Info.ProcessedPath)
% disp('===>> POSTPROCESSING FINISHED!')

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
           -Removed Limit of interpolation curve (AGr?ssl, Matthias Dieringer)
            that cut off high B1 values.
           -Automated the offset of the interpolation curve. Factors of
            360? are added automatically now (Matthias Dieringer)
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
            -optional input argument enSaving: 1 -> data and images are beeing saved
                                               0 -> no saving
            (Olli Kraus)
2014-05-06 -no Sensitivity.tif will be saved
            -enSaving (logical) is now replaced by SavingString (string). Empty string means no saving
            - allData.mat now is renamed to SavingString_MIDxxx.mat
            (Olli Kraus)
2016_06-20 -some errors in B1-calculation were fixed
           -external function calculates the exact RF pulse shape,
            which is needed to calculate the B1map from the FAmap
            -statistics section enabled
           - output is now a structure: Maps
2016_07_05 - modified to acurately calculate FA from 0?<=FA<=180? => extended DAM
            - before it was only 0?<=FA<=90?. For more information look into:
                An analysis of the accuracy of magnetic resonance flip
                angle measurement methods
                Glen R Morrell1 and Matthias C Schabel
                doi:10.1088/0031-9155/55/20/008
2016_07_21 - function can now readin rfpulsefile.txt which must be simulated with POET and saved as .txt file
2016_08_23 - now ether simulated, calculaed (assume: hanning-filtered sinc) or measured RF pulse shape are used for B1_abs calculation.

%}
end
