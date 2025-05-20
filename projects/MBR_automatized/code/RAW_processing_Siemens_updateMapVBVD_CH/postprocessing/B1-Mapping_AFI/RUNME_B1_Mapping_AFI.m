%% B1 Mapping based on Actual Flip Angle Imaging (AFI)
%--------------------------------------------------------------------------
% Input:
% 1. AFI-Folder
% 2. Reference Voltage
% Output: Images B1
%--------------------------------------------------------------------------
% Info:
% CONFIDENTIAL!!! NOT OUTSIDE B.U.F.F.:
% Flip Angle Aquirement is performed by Sebastian Schmitter (PTB)
%--------------------------------------------------------------------------
% Info:
% Author: Thomas Eigentler 
% Date: 2018-12-07
%--------------------------------------------------------------------------
% Trackchanges:
% 2018-12-07: Final Version 1.0 (Thomas Eigentler)

clear all;
close all

%% Reconstruct the Flipangle-Images (DICOM) @Sebastian Schmitter PTB

AFI = AFI_read_TXfct()                                                      %FlipAngle Map aquired by dividing 2 different contrast immages of the same slice

config.SetFA = AFI.mrprot.adFlipAngleDegree;                                %[?] set FA
[Export.Path, Export.FileName, Export.FileExtention] = fileparts(AFI.sDcmInfo.Filename);


%% Settings

config.RefVoltage = AFI.mrprot.sTXSPEC.asNucleusInfo.flReferenceAmplitude;  %[V]
config.PulsDuration = 1;                                                    %[ms] default setting of the Sequence
config.Z0 = 50;                                                             %[Ohm]
config.gamma = 2.6752219E8;                                                 %[rad/s/T]
config.TargetFAdeg = 180;                                                   %[?] default setting of the Sequence
config.TargetFArad = pi/180*config.TargetFAdeg;

%% Calculate the B1+ Map

B1.AngleDeg = AFI.alpha;                                                    %[?]
B1.RefVoltMap = config.SetFA/B1.AngleDeg.*config.RefVoltage                 %[V]
B1.Ref = config.TargetFArad/(config.gamma*config.PulsDuration*1e-3);        %[T] B1 for 180? FA @ Puls.Duration = 1ms
B1.Map = B1.Ref./sqrt(B1.RefVoltMap.^2./config.Z0).*sqrt(1000);             %[T/sqrt(kW)]
B1.Map = B1.Map.*1e6;                                                       %[uT/sqrt(kW)]     

%% Plot

SliceSelection = 12
FontSize = 20;
FigureSize = 800;

close all;
WriteImage = true;

%write B1 Map
ExportName = [Export.Path '\B1_Data.mat'];
save(ExportName,'B1')

%Plot Magnitude Image
i = figure(1)
    set(gcf,'Position',[0 0 FigureSize FigureSize])
subplot(2,1,1)
    imagesc(abs(AFI.DCM_TR1(:,:,SliceSelection)));
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial'); 
    colormap gray;
    c = colorbar;
    txt = ['Magnitude Image TR1 = ' num2str(AFI.TR1/1000) ' ms'];
    t1 = title(txt,'FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');
    axis image;
    shg;   
subplot(2,1,2)
    imagesc(abs(AFI.DCM_TR2(:,:,SliceSelection)));
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial'); 
    colormap gray;
    c = colorbar;
    txt = ['Magnitude Image TR2 = ' num2str(AFI.TR2/1000) ' ms'];
    t1 = title(txt,'FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');
    axis image;
    shg;
    if WriteImage
        ExportName = [Export.Path '\MagnitudeImage_Slice' num2str(SliceSelection)];
        saveas(i,ExportName,'fig')
        saveas(i,ExportName,'pdf')
        saveas(i,ExportName,'jpg')
    end

%Plot FA map
i = figure(2)
    set(gcf,'Position',[0 0 FigureSize FigureSize])
    imagesc(abs(B1.AngleDeg(:,:,SliceSelection))); 
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial'); 
    colormap jet;
    c = colorbar;
    caxis([0 90]);
    c.Label.String = ('$\mathrm{[^\circ]}$');
    set(c.Label,'Interpreter','latex');
    set(c.Label,'FontName','Arial');
    set(c.Label,'FontSize',FontSize);
    t1 = title('Actual Flip Angle Imaging (AFI): $\mathrm{FA-Map}$','FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');
    axis image;
    shg;
    if WriteImage
        ExportName = [Export.Path '\FA-Map_Slice' num2str(SliceSelection)];
        saveas(i,ExportName,'fig')
        saveas(i,ExportName,'pdf')
        saveas(i,ExportName,'jpg')
    end

%Plot B1+ Map
i = figure(3)
    set(gcf,'Position',[0 0 FigureSize FigureSize])
    imagesc(abs(B1.Map(:,:,SliceSelection))); 
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial'); 
    colormap jet;
    c = colorbar;
    caxis([0 50]);
    c.Label.String = ('$\mathrm{\mu T / \sqrt{kW}}$');
    set(c.Label,'Interpreter','latex');
    set(c.Label,'FontName','Arial');
    set(c.Label,'FontSize',FontSize);
    t1 = title('Actual Flip Angle Imaging (AFI): $\mathrm{B _{1} ^{+}}$','FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');
    axis image;
    shg;
    if WriteImage
        ExportName = [Export.Path '\B1-Map_Slice' num2str(SliceSelection)];
        saveas(i,ExportName,'fig')
        saveas(i,ExportName,'pdf')
        saveas(i,ExportName,'jpg')
    end

    
    
    
    
    
 
    
     
    
