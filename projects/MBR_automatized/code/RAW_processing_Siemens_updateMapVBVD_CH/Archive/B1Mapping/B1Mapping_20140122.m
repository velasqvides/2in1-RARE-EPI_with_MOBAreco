function [B1Map,B1Map_abs,B0Map]=B1Mapping_20140122(combinedImage,newDataDims,MrProt,config)
% This function calculates B0, normalized B1+ maps, and absolute B1+ maps for
% various phase based methods such as Bloch-Siegert-Shift, adiabatic Bloch-
% Siegert Shift, PhiFA-Cup and some more
%
% INPUT: combinedImage = complex images with four echoes
%        newDataDims = cell array with description of the data dimensions
%        MrProt = Struct with MR protocol parameters
%        config = struct with information such as filename and pathname
%
% OUTPUT: B1Map = B1+ map normalized to the input B1+
%         B1Map_abs = absolute B1+ map in [uT/sqrt(kw)]
%         B0Map = B0 map in [Hz]
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

% first we need to transform the data to double, since all other function
% use doubles. We would need to rewrite the whole program, especially the
% sensitivity calculation, otherwise
combinedImage=cast(combinedImage,'double');

Temp.Path=config.pathname;

Temp.MID=regexpi(config.filename,'_','split');
Temp.MID=Temp.MID{2};
Temp.Path=strcat(Temp.Path,'Processed_',Temp.MID,'_',date);

mkdir(Temp.Path);
Temp.Path1=fullfile(Temp.Path,'B0map');
mkdir(Temp.Path1);
Temp.Path2=fullfile(Temp.Path,'B1map');
mkdir(Temp.Path2);
Temp.Path3=fullfile(Temp.Path,'Data');
mkdir(Temp.Path3);
Temp.Path4=fullfile(Temp.Path,'Pictures');
mkdir(Temp.Path4);

% for the B1+ Mapping routine a certain order of the indices is more convenient
CellOrder={'Col','Lin','Eco','Par','Sli'};

% delete 'Channels' from the dimensions because the channels were already
% combined
newDataDims(strcmp(newDataDims,'Cha'))=[];

%get the permutation indices to get to the new order
for n=1:length(CellOrder)
    CellPos=find(~cellfun(@isempty,strfind(newDataDims,CellOrder{n})));
    if ~isempty(CellPos)
        permOrder(n)=CellPos;
    end
end

%permute the data according to the new order
CmplxImage=permute(combinedImage,permOrder);

% extract dimensions into variables
nPar=size(CmplxImage,4); % number of 3D partitions/slices

% get meta data
[Pulse,Sequence]=getRAWmetadata(MrProt);

% calculate phase differences for B0 map
imPhaDiffB0=squeeze(mod(angle(CmplxImage(:,:,2,:))-angle(CmplxImage(:,:,1,:))+pi,2*pi)-pi);

% calculate phase differences for B1 map
% alternatively, use echo 2 and 4 if 1 and 3 have artifacts
imPhaDiffB1=squeeze(mod(angle(CmplxImage(:,:,1,:))-angle(CmplxImage(:,:,3,:))+pi,2*pi)-pi);

% calculate magnitude image and scale to 4095 (this is arbitrary but
% similar to DICOM; remember this is single precision float!)
Magnitude1=squeeze(abs(CmplxImage(:,:,1,:)));
Magnitude1=Magnitude1/max(max(max(Magnitude1)))*4095;
disp('===>> Reconstruction done!');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    %%%
%%%  PHASE UNWRAPPING  %%%
%%%                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for cPar=1:nPar % display B0 phase differences to check for phase wraps
    Image_plot3D(1,nPar,cPar,0,1,'gray',imPhaDiffB0(:,:,cPar),0)
end
for cPar=1:nPar % display B1 phase differences to check for phase wraps
    Image_plot3D(2,nPar,cPar,0,1,'gray',imPhaDiffB1(:,:,cPar),0)
end

% Show magnitude and phase (does everything look good?) (optional)
% does your phase data need unwrapping?
disp('Please check B0 and B1 phase differences for phase wraps!');
inpt=input('Does your data need phase unwrapping [1] or not [Enter]?');
close all
if (~isempty(inpt))
    for cPar=1:nPar
        subplot(1,2,1)
        imagesc(imPhaDiffB0(:,:,cPar))
        title('Phase difference used for B_0-Map')
        axis image
        subplot(1,2,2)
        imagesc(imPhaDiffB1(:,:,cPar))
        title('Phase difference used for B_1-Map')
        axis image
        colormap(gray(256))
        inpt=input('Are there relevant phase wraps in B0 [1], B1 [2], both [3] or none [Enter]?');
        
        % perform unwrapping
        if (~isempty(inpt))
            cmplx1=abs(CmplxImage(:,:,1,cPar)).*exp(1i.*imPhaDiffB0(:,:,cPar));
            cmplx2=abs(CmplxImage(:,:,1,cPar)).*exp(1i.*imPhaDiffB1(:,:,cPar));
            
            while 1
                if (inpt==1 || inpt==3)
                    close all
                    imagesc(angle(cmplx1))
                    colorbar
                    axis image
                    title('Pick a point near (real) zero phase and with good signal!')
                    [xpoint,ypoint] = ginput(1);
                    disp('===>> Unwrapping, please wait...')
                    imPhaDiffB0(:,:,cPar)=Unwrap2D_QualityGuided_inputXY(cmplx1, 0.01 ,0,xpoint,ypoint);
                    disp('===>> Unwrapping done!')
                end
                
                if (inpt==2 || inpt==3)
                    close all
                    imagesc(angle(cmplx2))
                    colorbar
                    axis image
                    title('Pick a point near (real) zero phase and with good signal!')
                    [xpoint,ypoint] = ginput(1);
                    disp('===>> Unwrapping, please wait...')
                    imPhaDiffB1(:,:,cPar)=Unwrap2D_QualityGuided_inputXY(cmplx2, 0.01 ,0,xpoint,ypoint);
                end
                
                subplot(1,2,1)
                imagesc(imPhaDiffB0(:,:,cPar))
                colormap gray
                title('Phase difference used for B_0-Map')
                axis image
                colorbar
                subplot(1,2,2)
                imagesc(imPhaDiffB1(:,:,cPar))
                title('Phase difference used for B_1-Map')
                axis image
                colorbar
                inpt2=input('DONE! Happy with the result [any key] or not [1]?');
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
clear B0 B1
B0.Sub=imPhaDiffB0*180/pi;
B1.Sub=imPhaDiffB1*180/pi;

% cleanup workspace
clearvars -except Temp B0 B1 Magnitude1 Pulse Sequence nPar %clean up workspace
close all

%Display meta data
Pulse
Sequence


%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%           %%%
%%%  MASKING  %%%
%%%           %%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

figure
[Mask.Mask,Mask.Perim,Mask.Image]=masking_best_EM3D(Magnitude1);
close(gcf)
for cPar=1:nPar
    Image_plot3D(1,nPar,cPar,0,1,'gray',Mask.Mask(:,:,cPar),0)
end
set(gcf, 'name', 'Mask')
Image_save(fullfile(Temp.Path4,'Mask'),'-dtiff')
close(gcf)

% apply mask
B0.Sub=B0.Sub.*Mask.Mask;
B1.Sub=B1.Sub.*Mask.Mask;


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  B0 CALCULATIONS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate B0 Map
disp('===>> Calculating B0 maps...');
B0.Map=B0.Sub/(360*Sequence.DeltaTE);
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-180,180,'jet',B0.Map(:,:,cPar),1)
end
set(gcf, 'name', 'B0 map')

%cut off outliers
Temp.B0cut=input(['Maximum B0 is ',num2str(round(max(abs(B0.Map(:))))), 'Hz, choose cut off value ([Enter] to use max)> ']);
Temp.B0cut(isempty(Temp.B0cut))=max(abs(B0.Map(:)));
B0.Map(abs(B0.Map)>Temp.B0cut)=NaN;

% save image and data
close gcf
for cPar=1:nPar
    Image_plot3D(2,nPar,cPar,-max(abs(B0.Map(:))),max(abs(B0.Map(:))),'jet',B0.Map(:,:,cPar),1)
end
set(gcf, 'name', 'B0 map [Hz]')
Image_convert(2,Temp.Path1,'B0map')
save(fullfile(Temp.Path3,'Maps'),'Mask','B0','-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  B1 CALCULATIONS  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===>> Calculating B1+ maps...');

%%%%%%%%%%%%%%%%%%%%%
%%%  Sensitivity  %%%
%%%%%%%%%%%%%%%%%%%%%

%round up B0 distribution and add +20Hz
Data.DistrB0=ceil(max(abs(B0.Map(:)))/10)*10;

%number of spins adapted to data
Data.SpinsB0=Data.DistrB0/10;% this is good for Bloch Siegert
Data.B1Max=4;

%number of spins adapted
Data.SpinsB1=30*Data.B1Max;
Data

% non-mex version is not supported anymore in this version
Pulse.Shift=0;
% calculate sensitivity curve
[Curve]=Sensitivity_build_mex(Pulse,Data);
save(fullfile(Temp.Path3,'Sensitivity'),'Curve','Data','-v7.3')

figure(3)
surf(Curve.B0mat,Curve.B1mat,Curve.DeltaPhi)
axis square
set(gca,'xLim',[-Data.DistrB0-50 Data.DistrB0+50],'yLim',[0 max(Curve.B1mat(:))],'zLim',[min(Curve.DeltaPhi(:))-10 max(Curve.DeltaPhi(:))+10],'View',[-35 25],'FontSize',16)
zlabel('Phase accrual \psi [deg]'); ylabel('B_1 [normalized]'); xlabel('B_0 offset [Hz]')
set(gcf, 'name', 'B1+ sensitivity map: Phase accrual')
Image_save(fullfile(Temp.Path4,'Sensitivity'),'-dtiff')


%%%%%%%%%%%%%%%%%%%%%%%
%%%  Interpolation  %%%
%%%%%%%%%%%%%%%%%%%%%%%

B0.Map(isnan(B0.Map))=0; %cut out possible NaNs from phase unwrapping and thresholding
for cPar=1:nPar
    B1.Map(:,:,cPar)=Mask.Mask(:,:,cPar).*Interpolation_B1B0(B0.Map(:,:,cPar),B1.Sub(:,:,cPar),Curve.DeltaPhi,Curve.B0vec,Curve.B1vec);
end
Mask.Effective=Mask.Mask;
Mask.Effective(B1.Map==0)=0;


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%              %%%
%%%  STATISTICS  %%%
%%%              %%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

disp('===>> Calculating statistics...');
Temp.Bins=(0.005:0.01:Data.B1Max);

Slice.B1avg(:,1)=sum(sum(B1.Map,1),2)./sum(sum(Mask.Effective,1),2);
Slice.B1avg(isnan(Slice.B1avg))=0;
Slice.B1avg_avg=sum(Slice.B1avg)/sum(Slice.B1avg~=0);
Slice.B1avg_std=sqrt(sum((Slice.B1avg-Slice.B1avg_avg).^2)/sum(Slice.B1avg~=0));

Slice.B1std=zeros(nPar,1);
for cPar=1:nPar
    Slice.B1std(cPar,1)=sqrt(sum(sum((B1.Map(:,:,cPar)-Mask.Effective(:,:,cPar)*Slice.B1avg(cPar,1)).^2,1),2)/sum(sum(Mask.Effective(:,:,cPar),1),2));
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

figure(4)
bar(Temp.Bins,Total.B1dist);
axis square
set(gcf, 'name', 'B1+ histogram')
set(gca,'XLim',[0 max(Temp.Bins)])
title('B_1^+ distribution')
ylabel('Normalized Distribution')
xlabel(['AVG = ' num2str(Total.B1dist_avg) '   STD = ' num2str(Total.B1dist_std)])
Image_save(fullfile(Temp.Path4,'B1dist'),'-dtiff')


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    %%%
%%%  DISPLAY B1+ MAPS  %%%
%%%                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% normalized B1+ MAP %
%%%%%%%%%%%%%%%%%%%%%%
for cPar=1:nPar
    Image_plot3D(5,nPar,cPar,0,2,'hot',B1.Map(:,:,cPar),1)
end
set(gcf, 'name', 'normalized B1+ [%/100]')
Image_convert(5,Temp.Path2,'B1map')
B1Map=B1.Map;
B0Map=B0.Map;
disp(['Mean normalized B1+ = ',num2str(100*mean2(B1Map(B1Map>0))),' +- ',num2str(100*std2(B1Map(B1Map>0))),' %'])
save(fullfile(Temp.Path,'B0B1Map_from_RAW.mat'),'B0Map','B1Map','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% absolute B1+ map in uT/sqrt(kW) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1Map_abs = B1.Map.*Pulse.NominalB1./Sequence.ReqVolt.*sqrt(50).*sqrt(1000);
for cPar=1:nPar
    Image_plot3D(6,nPar,cPar,0,10,'hot',B1Map_abs(:,:,cPar),1)
end
set(gcf, 'name', 'Absolute B1+ [uT/sqrt(kW)]')
Image_convert(6,Temp.Path2,'B1map_abs')
save(fullfile(Temp.Path,'B0B1Map_from_RAW.mat'),'B1Map_abs','-v7.3','-append')
save(fullfile(Temp.Path3,'Maps'),'Mask','B0','B1','B1Map_abs','-v7.3')
disp('===>> EVALUATION FINISHED!')

% changelog:
%
% 2013-05-16 -Adapted from Complete_F_B1new.m to work with Siemens RAW data
%             (Matthias Dieringer)
% 2013-05-22 -added averages
%             (Matthias Dieringer)
% 2013-05-21 -Automatic image registration implemented
%             (Matthias Dieringer)
% 2013-05-29 -changed order of cells to calculate everything with one click
%            -Flaw with registration fixed
%             (Matthias Dieringer)
% 2013-05-30 -added phase unwrapping
%            -Removed Limit of interpolation curve (AGrässl, Matthias Dieringer)
%             that cut off high B1 values.
%            -Automated the offset of the interpolation curve. Factors of
%             360° are added automatically now (Matthias Dieringer)
%            -dramatic speed increase for calculation of the interpolation
%             curve. If you don't have the correct MEX-file or can't
%             compile it, the calculation will switch back to the slower
%             Matlab version and ask for desired B1 scale factor
%             (Matthias Dieringer)
% 2013-06-06 -Speed-up and simplification in the Interpolation_B0B1
%             function. There were sometimes problems with the lookup
%             table (leading to extremely high B1 in some regions),
%             which should be fixed now.
%             (Matthias Dieringer)
% 2013-06-13 -changed raw reading routine
%            -extended raw data filtering options
%            -code is faster and simplified
%             (Matthias Dieringer)
% 2013-06-13 -adapted code to handle 1 receiver only as well
%             (Matthias Dieringer)
% 2013-06-17 -RF pulse and sequence info is fetched from RAW data
%             header (Matthias Dieringer)
% 2013-06-17 -removed some unnecessary code and graphs
%            -changed max B1 sensitivity (curve) to 4 with 120 samples
%            -changed how the B0 Map is cut off (was with sign switch
%             before, now outliers are set to 0)
%             (Matthias Dieringer)
% 2013-07-02 -speed improvements and simplifications in sensitivity_build
%             and pulse_build
%            -non-MEX version of sensitivity build is not supported anymore
%             in this version
%            -Adiabatic Bloch Siegert pulse implemented (Khaligi et al, MRM
%             2013)
%            -Gaussian RF pulse was removed, no application for this
%            (Matthias Dieringer)
%
%            THEREFORE YOU NEED TO USE THE NEW VERSION OF THE SEQUENCE
%            (20130702 and later)
%
% 2013-08-29 -there was a bug in the Pulse_build.m so that the scaling was
%             not correct (only the new version 2013-07-02 is affected).
%             This could lead to an B1 error of up to 5%
%             (Matthias Dieringer)
% 2013-09-03 -change save(...) to save(..,'-v7.3') to allow huge files and
%             to allow loading of single variables out of a saved array
%            -change strcat to fullfile in order to make it working for
%             Linux as well
%            -B0 distribution to simulate is rounded up instead of being
%             rounded down (caused NaNs in the interpolation sonmetimes)
%             (Olli Kraus)
% 2013-09-06 -replaced single precision variables by double precision
%             variables to make this script compatible for use with Linux
%            -Linux 64bit mex bloch simulation file added
% 2013-09-25 -Added absolute Map in ut/(sqrt/kW). Reference point for the
%             RMS input power is the siemens coil plug(A Graessl)
% 2013-09-27 -meas ID is added to the name of the subfolder "Processed" and
%             maps are stored in the "processed" subfolder.
%            -after unwrapping, "any key" is used to proceed instead of "1"
%            -new function "Display_Image" to display images
%            -cleaned up code
%            (Matthias Dieringer)
% 2014-01-13 -3D and multislice is enabled again
%            -3D methods using composite RF hard pulses are facilitated now
%            -RAW data import and handling adapted
%            -MR Sequence was adapted as well, please use the new version
%            (although older RAW files should still work)
%            -cleaned up some of the code (also external functions)
end