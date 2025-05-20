function TemperaturePlotter(TempMap, ExportPath, ExportFileName, MaxTemp, FontSize, FigureSize)

load('ColormapCST.mat')

i = figure(1)
set(gcf,'Position',[0 0 FigureSize FigureSize])
imagesc(TempMap,[0 MaxTemp]);
set(gca,'FontSize',FontSize); 
set(gca,'FontName','Arial'); 
%colormap(ColormapCST);
colormap(parula);
c = colorbar;
t1 = title(ExportFileName,'FontSize', FontSize);
y1 = ylabel(colorbar,'\DeltaT [\?C]','FontSize', FontSize);
axis image;
shg;   

ExportName = [ExportPath, ExportFileName];
saveas(i,ExportName,'fig');
saveas(i,ExportName,'pdf');
saveas(i,ExportName,'jpg');

end

