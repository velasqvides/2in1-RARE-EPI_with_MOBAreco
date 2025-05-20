function [B1Map_rel,B1Map_uT,B0Map_Hz,detailedData]=B1Mapping_20140221(combinedImage,newDataDims,MrProt,config)
% This function calculates B0, normalized B1+ maps, and absolute B1+ maps for
% various phase based methods such as Bloch-Siegert-Shift, adiabatic Bloch-
% Siegert Shift, PhiFA-Cup and some more
%
% INPUT: combinedImage = complex images with four echoes
%        newDataDims = cell array with description of the data dimensions
%        MrProt = Struct with MR protocol parameters
%        config = struct with information such as filename and pathname
%
% OUTPUT: B1Map_rel = relative B1 map (wherever B1Map_rel=100% => actuaal FA = nominal FA of protocol)
%         B1Map_uT = absolute B1+ map in [uT/sqrt(kW)]
%         B0Map = B0 map in [Hz]
%         detailedData = structure of other relevant data
%
% changelog: see end of file
%
% 2014 Matthias Dieringer
% matthias.dieringer@charite.de

if(nargin~=4)
    disp('Wrong number of arguments!')
    help B1Mapping
    return
end

%% first we need to transform the data to double, since all other function
% use doubles. We would need to rewrite the whole program, especially the
% sensitivity calculation, otherwise
combinedImage=cast(combinedImage,'double');

Info.MID=regexpi(config.filename,'_','split');
Info.MID=Info.MID{2};
Info.Subfolder=strcat('Processed_',Info.MID,'_',date);
Info.Path=fullfile(config.pathname,Info.Subfolder);
mkdir(Info.Path);

% Temp.Path1=fullfile(Temp.Path,'B0map');
% mkdir(Temp.Path1);
% Temp.Path2=fullfile(Temp.Path,'B1map');
% mkdir(Temp.Path2);
% Temp.Path3=fullfile(Temp.Path,'Data');
% mkdir(Temp.Path3);
% Temp.Path4=fullfile(Temp.Path,'Pictures');
% mkdir(Temp.Path4);

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

% get meta data
[Pulse,Sequence]=getRAWmetadata(MrProt);

% calculate phase differences for B0 map
%imPhaDiffB0=squeeze(mod(angle(CmplxImage(:,:,2,:))-angle(CmplxImage(:,:,1,:))+pi,2*pi)-pi);
imPhaDiffB0=squeeze(phaseDiff(CmplxImage(:,:,2,:),CmplxImage(:,:,1,:)));

% calculate phase differences for B1 map
% alternatively, use echo 2 and 4 if 1 and 3 have artifacts
%imPhaDiffB1=squeeze(mod(angle(CmplxImage(:,:,1,:))-angle(CmplxImage(:,:,3,:))+pi,2*pi)-pi);
imPhaDiffB1=squeeze(phaseDiff(CmplxImage(:,:,1,:),CmplxImage(:,:,3,:)));

% calculate magnitude image and scale to 4095 (this is arbitrary but
Magnitude1=squeeze(abs(CmplxImage(:,:,1,:)));
%max(Magnitude1(:))
%Magnitude1=Magnitude1/max(Magnitude1(:));%normalized to 1
%disp('===>> Reconstruction done!');

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
[Masks.Mask,Masks.Image,Info.relativeMagnitudeThreshold]=masking_by_thresholding_3D(Magnitude1);

close(gcf)
for cPar=1:nPar
    Image_plot3D(1,nPar,cPar,0,1,'jet',Masks.Mask(:,:,cPar),0,1)
end
set(gcf, 'name', 'ThresholdMask')
Image_save(fullfile(Info.Path,'ThresholdMask'),'-dtiff')
%close(gcf)
%% optional: draw ROI for every slice

Info.enROI=input('Do you want to draw a ROI [1] or not [Enter]? ');
%close all
if  (~isempty(Info.enROI))
    for cPar=1:nPar
        figure(11)
        imagesc(Masks.Image(:,:,cPar));
        figID=gca;   
        ROI=imfreehand(figID);
        %ROI=imellipse(figID);
        %ROI=imrect(figID);
        %wait(ROI);
        Masks.ROI(:,:,cPar)=createMask(ROI);
        Image_plot3D(16,nPar,cPar,0,1,'jet',Masks.ROI(:,:,cPar).* Masks.Image(:,:,cPar),0,1)    
    end
    set(gcf, 'name', 'ROIMask')
    Image_save(fullfile(Info.Path,'ROIMask'),'-dtiff')
else
    Info.enROI=0;
    disp('no ROI drawn')
end
%close all

%% %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    %%%
%%%  PHASE UNWRAPPING  %%%
%%%                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for cPar=1:nPar % display B0 phase differences to check for phase wraps
    data2plot=imPhaDiffB0(:,:,cPar).*Masks.Mask(:,:,cPar);
    Min=min(data2plot(:));
    Max=max(data2plot(:));
    Image_plot3D(21,nPar,cPar,Min,Max,'jet',data2plot,0,1)
end
for cPar=1:nPar % display B1 phase differences to check for phase wraps
    data2plot=imPhaDiffB1(:,:,cPar).*Masks.Mask(:,:,cPar);
    Min=min(data2plot(:));
    Max=max(data2plot(:));
    Image_plot3D(22,nPar,cPar,Min,Max,'jet',data2plot,0,1)
end

% Show magnitude and phase (does everything look good?) (optional)
% does your phase data need unwrapping?
disp('Please check B0 and B1 phase differences for phase wraps!');
Info.enUnwrapping=input('Does your data need phase unwrapping [1] or not [Enter]? ');
close all
if (~isempty(Info.enUnwrapping))
    for cPar=1:nPar
        subplot(1,2,1)
        imagesc(imPhaDiffB0(:,:,cPar).*Masks.Mask(:,:,cPar))
        set(gca,'Title',text('String',['Phase difference used for B_0-Map. slice: ', num2str(cPar)]));
        axis image
        subplot(1,2,2)
        imagesc(imPhaDiffB1(:,:,cPar).*Masks.Mask(:,:,cPar))
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
                    set(gca,'Title',text('String','Pick a point near (real) zero phase and with good signal!'))
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
                    set(gca,'Title',text('String','Pick a point near (real) zero phase and with good signal!'))
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
                %inpt2=input('DONE! Happy with the result [any key] or not [1]? ');
                inpt2=[];
                inpt2(isempty(inpt2))=0;
                if inpt2~=1
                    break
                end
            end
        end
    end
end

% prepare for export [deg]
close all
% cleanup workspace
clearvars -except Info Masks Magnitude1 Pulse Sequence nPar imPhaDiffB0 imPhaDiffB1 MrProt%clean up workspace

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
% B0.Sub=imPhaDiffB0*180/pi .*Masks.Mask;
% B0.Map=B0.Sub/(360*Sequence.DeltaTE);

B0.Map=imPhaDiffB0/2/pi/Sequence.DeltaTE .*Masks.Mask;
Info.MaxDeltaB0=max(abs(B0.Map(:)));
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-Info.MaxDeltaB0,Info.MaxDeltaB0,'jet',B0.Map(:,:,cPar),1,1)
end
set(gcf, 'name', 'B0 map')

%cut off outliers
disp(['Maximum B0 is ',num2str(round(Info.MaxDeltaB0)),' Hz']);
Info.B0cut=input('choose cut off value ([Enter] to use max)> ');
Info.B0cut(isempty(Info.B0cut))=Info.MaxDeltaB0;
B0.outliers=logical(abs(B0.Map)>Info.B0cut);
B0.noOutliers=sum(B0.outliers(:));
B0.Map(abs(B0.Map)>Info.B0cut)=0;

%% save image and data
close gcf
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-Info.B0cut,Info.B0cut,'jet',B0.Map(:,:,cPar).*(1-B0.outliers(:,:,cPar)),1,1)
end
set(gcf, 'name', 'B0 map [Hz]')
Image_save(fullfile(Info.Path,'B0Map'),'-dtiff')

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  B1 CALCULATIONS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===>> Calculating B1+ maps...');
%B1.Sub=imPhaDiffB1*180/pi .*Masks.Mask;

%%%%%%%%%%%%%%%%%%%%%
%%%  Sensitivity  %%%
%%%%%%%%%%%%%%%%%%%%%

%round up B0 distribution and add +20Hz
Data.DistrB0=ceil(max(abs(B0.Map(:)))/10)*10;

%number of spins adapted to data
Data.SpinsB0=Data.DistrB0/10;% this is good for Bloch Siegert
Data.B1Max=4; % 400% is a realistic upper limit. more is unrealistic
Data.SpinsB1=30*Data.B1Max;%number of spins adapted

% non-mex version is not supported anymore in this version
Pulse.Shift=0;
% calculate sensitivity curve
%[SimulatedSensitity]=Sensitivity_build_mex(Pulse,Data);
SimulatedSensitivity=Sensitivity_build_mex(Pulse,Data);

SimulatedSensitivity.B1Spins=Data.SpinsB1;
SimulatedSensitivity.B0Spins=Data.SpinsB0;
SimulatedSensitivity.B0Range=[-Data.DistrB0 +Data.DistrB0];
SimulatedSensitivity.B1Range=[0 Data.B1Max];
SimulatedSensitivity.DeltaPhiRange=[0 max(SimulatedSensitivity.DeltaPhi(:))];

%save(fullfile(Info.Path,'Sensitivity'),'SimulatedSensitivity','-v7.3')

figure(3)
surf(SimulatedSensitivity.B0mat,SimulatedSensitivity.B1mat,SimulatedSensitivity.DeltaPhi)
axis square
set(gca,'xLim',[-Data.DistrB0-50 Data.DistrB0+50],'yLim',[0 max(SimulatedSensitivity.B1mat(:))],'zLim',[min(SimulatedSensitivity.DeltaPhi(:))-10 max(SimulatedSensitivity.DeltaPhi(:))+10],'View',[-35 25],'FontSize',16)
zlabel('Phase accrual \psi [deg]'); ylabel('B_1 [normalized]'); xlabel('B_0 offset [Hz]')
set(gcf, 'name', 'B1+ sensitivity map: Phase accrual')
Image_save(fullfile(Info.Path,'Sensitivity'),'-dtiff')
clear Data

%% %%%%%%%%%%%%%%%%%%%%%
%%%  Interpolation  %%%
%%%%%%%%%%%%%%%%%%%%%%%

for cPar=1:nPar
    B1.Map_rel(:,:,cPar)=Masks.Mask(:,:,cPar).*Interpolation_B1B0(B0.Map(:,:,cPar),imPhaDiffB1(:,:,cPar)*180/pi .*Masks.Mask(:,:,cPar),SimulatedSensitivity.DeltaPhi,SimulatedSensitivity.B0vec,SimulatedSensitivity.B1vec);
end
B1.Map_rel(B0.outliers)=NaN;
B1.Singularities= ~isfinite(B1.Map_rel);
B1.noSingularities=sum(B1.Singularities(:));
Masks.Effective= ~B1.Singularities & Masks.Mask & ~B0.outliers;
if (Info.enROI)
    Masks.Effective_ROI = ~B1.Singularities & Masks.Mask & Masks.ROI & ~B0.outliers;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    %%%
%%%  DISPLAY B1+ MAPS  %%%
%%%                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized B1+ Map &  absolute B1+ map in uT/sqrt(kW) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1.rel2absFactor = Pulse.NominalB1 / Sequence.ReqVolt * sqrt(50) * sqrt(1000);
B1.Map_uT = B1.Map_rel * B1.rel2absFactor;%unit: uT/sqrt(kW)
for cPar=1:nPar
    Image_plot3D(4,nPar,cPar,0,2,'hot',B1.Map_rel(:,:,cPar),1,1)
end
set(gcf, 'name', 'normalized B1+ [%/100]')
Image_save(fullfile(Info.Path,'B1Map_rel'),'-dtiff')

Info.MaxB1_uT=max(B1.Map_uT(:));
for cPar=1:nPar
    Image_plot3D(5,nPar,cPar,0,Info.MaxB1_uT,'hot',B1.Map_uT(:,:,cPar),1,1)
end
set(gcf, 'name', 'Absolute B1+ [uT/sqrt(kW)]')
Image_save(fullfile(Info.Path,'B1Map_uT'),'-dtiff')

fprintf('Mean normalized B1 = %1.2f +- %1.2f %%\n',100*mean2(B1.Map_rel(B1.Map_rel>0)),100*std2(B1.Map_rel(B1.Map_rel>0)) );

%% %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%              %%%
%%%  STATISTICS  %%%
%%%              %%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%{
% old statistic calculation
disp('===>> Calculating statistics...');
Temp.Bins=(0.005:0.01:Data.B1Max);

Slice.B1avg(:,1)=sum(sum(B1.Map,1),2)./sum(sum(Masks.Effective,1),2);
Slice.B1avg(isnan(Slice.B1avg))=0;
Slice.B1avg_avg=sum(Slice.B1avg)/sum(Slice.B1avg~=0);
Slice.B1avg_std=sqrt(sum((Slice.B1avg-Slice.B1avg_avg).^2)/sum(Slice.B1avg~=0));

Slice.B1std=zeros(nPar,1);
for cPar=1:nPar
    Slice.B1std(cPar,1)=sqrt(sum(sum((B1.Map(:,:,cPar)-Masks.Effective(:,:,cPar)*Slice.B1avg(cPar,1)).^2,1),2)/sum(sum(Masks.Effective(:,:,cPar),1),2));
end
Slice.B1std(isnan(Slice.B1std))=0;
Slice.B1std_avg=sum(Slice.B1std)/sum(Slice.B1std~=0);
Slice.B1std_std=sqrt(sum((Slice.B1std-Slice.B1std_avg).^2)/sum(Slice.B1std~=0));

Temp.B1vec=B1.Map(:);
Temp.B1vec(Temp.B1vec<min(Temp.Bins))=[];
Temp.B1vec(Temp.B1vec>max(Temp.Bins))=[];
Total.B1dist=hist(Temp.B1vec,Temp.Bins);
Total.B1dist=Total.B1dist/(sum(Total.B1dist));
Total.B1dist(isnan(Total.B1dist))=0;
Total.B1dist_avg=sum(Total.B1dist.*Temp.Bins);
Total.B1dist_std=sqrt(sum(Total.B1dist.*(Temp.Bins-Total.B1dist_avg).^2));

Slice.B1dist=zeros(nPar,length(Temp.Bins));
for cPar=1:nPar
    Temp.B1vec=reshape(B1.Map(:,:,cPar),size(B1.Map,1)*size(B1.Map,2),1);
    Temp.B1vec(Temp.B1vec<min(Temp.Bins))=[];
    Temp.B1vec(Temp.B1vec>max(Temp.Bins))=[];
    Slice.B1dist(cPar,:)=hist(Temp.B1vec,Temp.Bins);
    Slice.B1dist(cPar,:)=Slice.B1dist(cPar,:)/sum(Slice.B1dist(cPar,:));
    Slice.B1dist(isnan(Slice.B1dist))=0;
end
save(fullfile(Temp.Path3,'Statistics'),'Slice','Total','-v7.3')
%}
Statistics.B1Data       =B1.Map_uT(Masks.Effective);
Statistics.B1Max        =max(Statistics.B1Data);
Statistics.Bins         =0 :Statistics.B1Max/(300-1) : Statistics.B1Max;% 1x300 vector
Statistics.B1hist       =histc(Statistics.B1Data,Statistics.Bins);
Statistics.B1hist_norm  =Statistics.B1hist/length(Statistics.B1Data);
Statistics.B1_avg       =mean(Statistics.B1Data);
Statistics.B1_std       =std(Statistics.B1Data);
Statistics.B1_CV        =Statistics.B1_std/Statistics.B1_avg;

figure(6)
bar(Statistics.Bins(2:end),Statistics.B1hist_norm(2:end));
axis square
set(gcf, 'name', 'B1+ histogram')
set(gca,'XLim',[0 Statistics.B1Max])
xlabel('B_1^+ distribution [uT]')
ylabel('Normalized Distribution')
set(gca,'Title',text('String',['AVG = ', num2str(Statistics.B1_avg), ' uT , STD = ' num2str(Statistics.B1_std),' uT']))
Image_save(fullfile(Info.Path,'B1hist'),'-dtiff')
%% %%%%%%%%%%%%%%%%%%%%%
% STATISTICS INSIDE ROI%
%%%%%%%%%%%%%%%%%%%%%%%%
% take only data inside ROI & above threshold & no B1 singularity
if (Info.enROI)
    Statistics.ROI.B1Data       =B1.Map_uT(Masks.Effective_ROI); 
    Statistics.ROI.B1Max        =max(Statistics.ROI.B1Data);
    Statistics.ROI.Bins         =0 : Statistics.ROI.B1Max/(300-1) : Statistics.ROI.B1Max;% 1x300 vector
    %Statistics.ROI.Bins         =(0:0.01:SimulatedSensitivity.B1Range(2));
    Statistics.ROI.B1hist       =histc(Statistics.ROI.B1Data,Statistics.ROI.Bins);
    Statistics.ROI.B1hist_norm  =Statistics.ROI.B1hist/length(Statistics.ROI.B1Data);
    Statistics.ROI.B1_avg       =mean(Statistics.ROI.B1Data);
    Statistics.ROI.B1_std       =std(Statistics.ROI.B1Data);
    Statistics.ROI.B1_CV        =Statistics.ROI.B1_std/Statistics.ROI.B1_avg;
    
    figure(7)
    bar(Statistics.ROI.Bins(2:end),Statistics.ROI.B1hist_norm(2:end));
    axis square
    set(gcf, 'name', 'B1+ histogram inside ROI')
    set(gca,'XLim',[0 Statistics.ROI.B1Max])
    xlabel('B_1^+ distribution inside ROI [uT]')
    ylabel('Normalized Distribution')
    set(gca,'Title',text('String',['AVG = ', num2str(Statistics.ROI.B1_avg),' uT, STD = ' num2str(Statistics.ROI.B1_std),' uT']))
    Image_save(fullfile(Info.Path,'B1hist_ROI'),'-dtiff')
    
    string1=sprintf('(B1act/B1nom)_ROI = %1.2f +- %1.2f uT',Statistics.ROI.B1_avg,Statistics.ROI.B1_std);
    %string2=sprintf('\nOld Reference Voltage: %3.1f V ',Sequence.RefVolt);
    string3=sprintf('\nAdjust Reference Voltage to %3.1f V ',Sequence.RefVolt/Statistics.ROI.B1_avg * B1.rel2absFactor);
    disp([string1,string3])
end

%% cleanup workspace, prepare export and save data
clearvars -except B0 B1 Info Magnitude1 Masks Pulse Sequence SimulatedSensitivity Statistics MrProt%clean up workspace
save(fullfile(Info.Path,'allData'),'-v7.3')

B1Map_rel=B1.Map_rel;
B0Map_Hz=B0.Map;
B1Map_uT=B1.Map_uT;
save(fullfile(Info.Path,['B0B1Maps_',Info.MID]),'B0Map_Hz','B1Map_rel','B1Map_uT','-v7.3')

detailedData=struct('Info',Info,'MrProt',MrProt,'Pulse',Pulse,'Sequence',Sequence,...
    'SimulatedSensitity',SimulatedSensitivity,'Magnitude_Echo1',Magnitude1,...
    'Masks',Masks,'B0',B0,'B1',B1,'Statistics',Statistics);

disp('data and images saved to: ')
disp(Info.Path)
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
%}
end
