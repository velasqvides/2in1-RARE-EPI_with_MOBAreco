function Image_save(fname,format)
%IMAGE_SAVE Summary of this function goes here
%   Detailed explanation goes here

set(gcf,'PaperPositionMode','auto');
if strcmpi(format,'fig')
    saveas(gcf,fname,'fig')
else
    print(gcf,'-r300',format,'-cmyk','-opengl',fname);
    % system(['convert ' fname '.tif' ' ' fname '.jpg']);
end
end
