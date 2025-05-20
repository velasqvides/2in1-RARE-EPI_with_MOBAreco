%% B1 Mapping based on WIP543
%--------------------------------------------------------------------------
% Input: NON
% Output: Images B1
%--------------------------------------------------------------------------
% mrprot:
% Author: Thomas Eigentler 
% Date: 2018-10-05
%--------------------------------------------------------------------------
% Trackchanges:
% 2018-10-05: Final Version 1.0 (Thomas Eigentler)

clear all;
close all

%% Data Read

[Image.Data, Image.mrprot, Image.Pathname, Image.Filenames] = readDcmFolder()

[Export.Path, Export.FileName, Export.FileExtention] = fileparts(Image.mrprot.Filename);

%% Settings
% This sequence is only accurate with Rectangular Puls!

% Reading reference Data of the location
clear config;
config.pathname = 'P:\Measurement_7T\20181103_CardiacArraySGBT_B1Mapping\Cardiac_Array\MDC-0169_20181103_B1Mapping1\MDC-0169_20181103_B1Mapping1\1_TE^B1_Mapping_ID00004996';

config.RefVoltage = 100;                                                    %[V] VRMS for 180? FA @ Puls.Duration Reference Voltage
config.Z0 = 50;                                                             %[Ohm]
config.PulsDuration = 1;                                                    %[ms] default Setting for the sequence
config.gamma = 2.6752219E8;                                                 %[rad/s/T]
config.TargetFAdeg = 180;                                                   %[?] default setting of the Sequence
config.TargetFArad = pi/180*config.TargetFAdeg;

config.SetFA = Image.mrprot.FlipAngle;                                      %[?] set FA

%% Flipangle Map and B1+ field calculation

B1.AngleDeg = double(Image.Data)/10;                                        %[?]
B1.RefVoltMap = config.SetFA./B1.AngleDeg.*config.RefVoltage                %[V]
B1.Ref = config.TargetFArad/(config.gamma*config.PulsDuration*1e-3);        %[T] B1 for 180? FA @ Puls.Duration
B1.Map = B1.Ref./sqrt(B1.RefVoltMap.^2./config.Z0).*sqrt(1000);             %[T/sqrt(kW)]
B1.Map = B1.Map.*1e6;                                                       %[uT/sqrt(kW)]      

%% Plot 

%% Plot

FontSize = 20;
FigureSize = 800;

WriteImage = true;

%Plot Magnitude Image
i = figure(1)
    set(gcf,'Position',[0 0 FigureSize FigureSize])
    imagesc(abs(Image.Data));
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial'); 
    colormap gray;
    c = colorbar;
    t1 = title('Presaturation Based Mapping (WIP543): Magnitude Image','FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');
    axis image;
    shg;   
    if WriteImage
        ExportName = [Export.Path '\MagnitudeImage'];
        saveas(i,ExportName,'fig')
        saveas(i,ExportName,'pdf')
        saveas(i,ExportName,'jpg')
    end

%Plot FA map
i = figure(2)
    set(gcf,'Position',[0 0 FigureSize FigureSize])
    imagesc(abs(B1.AngleDeg)); 
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial');  
    colormap jet;
    c = colorbar;
    caxis([0 90]);
    c.Label.String = ('$\mathrm{[^\circ]}$');
    set(c.Label,'Interpreter','latex');
    set(c.Label,'FontName','Arial');
    set(c.Label,'FontSize',FontSize);
    t1 = title('Presaturation Based Mapping (WIP543): $\mathrm{FA-Map}$','FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');
    axis image;
    shg;
    if WriteImage
        ExportName = [Export.Path '\FA-Map'];
        saveas(i,ExportName,'fig')
        saveas(i,ExportName,'pdf')
        saveas(i,ExportName,'jpg')
    end

%Plot B1+ Map
i = figure(3)
    set(gcf,'Position',[0 0 FigureSize FigureSize])
    imagesc(abs(B1.Map)); 
    set(gca,'FontSize',FontSize); 
    set(gca,'FontName','Arial'); 
    colormap jet;
    c = colorbar;
    caxis([0 40]);
    c.Label.String = ('$\mathrm{\mu T / \sqrt{kW}}$');
    set(c.Label,'Interpreter','latex');
    set(c.Label,'FontName','Arial');
    set(c.Label,'FontSize',FontSize);
    t1 = title('Presaturation Based Mapping (WIP543): $\mathrm{B _{1} ^{+}}$','FontSize', FontSize);
    set(t1,'Interpreter','Latex')
    set(t1,'FontName','Arial');    
    axis image;
    shg;
    if WriteImage
        ExportName = [Export.Path '\B1-Map'];
        saveas(i,ExportName,'fig')
        saveas(i,ExportName,'pdf')
        saveas(i,ExportName,'jpg')
    end

