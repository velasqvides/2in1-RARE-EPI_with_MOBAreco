function Image_convert(FigID,path,fname)
%IMAGE_CONVERT Summary of this function goes here
%   Detailed explanation goes here

figure(FigID)
Image_save(strcat(path,'/',fname),'-dtiff')
system(['convert -loop 0 -delay 50 ' path '/Partition*.jpg' ' ' path '/' fname '.gif']);

end

