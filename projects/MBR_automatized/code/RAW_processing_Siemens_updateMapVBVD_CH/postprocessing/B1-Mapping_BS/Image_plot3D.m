function Image_plot3D(FigID,Partitions,i,MinVal,MaxVal,Colors,Map,Bar,enMaskedPlotting)
%IMAGE_PLOT3D Summary of this function goes here
%   Detailed explanation goes here
if nargin<9
    enMaskedPlotting=0;
else
    if ischar(Colors)
        Colors=eval(Colors); % Colors is now a M x 3 array that defines the color
    end
end

if(Partitions>1) %3D stack
    pixels_axes=120;
else %only one image
    pixels_axes=320;
end
pixels_bar=100;

size_ima1=ceil(Partitions/ceil(sqrt(Partitions)));
size_ima2=ceil(sqrt(Partitions));
index1=ceil(i/size_ima2);
index2=i-(index1-1)*size_ima2;

figure(FigID)
if i==1
    set(gcf,'Units','pixels','Position',[200 150 pixels_axes*size_ima2+pixels_bar*Bar pixels_axes*size_ima1]);
end
axes('Units','pixels','Position',[pixels_axes*(index2-1)+1 pixels_axes*(size_ima1-index1)+1 pixels_axes pixels_axes]);

% here is the actual plot-command
if ~enMaskedPlotting
    imagesc(Map,[MinVal MaxVal]);
else
    Map(isnan(Map))=0;
    plotDataRGB=RGBsetZerosAsWhite(Map,Colors,[MinVal,MaxVal]);
    imagesc(plotDataRGB,[MinVal MaxVal]);
    %imagescwithnan(Map,Colors,[1 1 1]);
end

if Bar==1&&i==Partitions
%    colorbar('Units','pixels','FontSize',16,'OuterPosition',[pixels_axes*size_ima2+pixels_bar/5 1 pixels_bar-pixels_bar/5 pixels_axes*size_ima1]);
% TE ge??ndert
    set(gca,'Position',[pixels_axes*(index2-1)+1 pixels_axes*(size_ima1-index1)+1 pixels_axes pixels_axes]);
end
set(gca,'xTickLabel',[],'xTick',[],'yTickLabel',[],'yTick',[]);
colormap(Colors);

end
