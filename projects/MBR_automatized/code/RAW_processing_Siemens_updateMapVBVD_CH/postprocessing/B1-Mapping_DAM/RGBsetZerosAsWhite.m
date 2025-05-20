function outRGB=RGBsetZerosAsWhite(inputdata,colorScheme,Clims)


imageRGB=mapcolor(inputdata,colorScheme,Clims);
mask=logical(inputdata); % sets only values which are exactly 0 to 0; rest is 1

[noX noY]=size(inputdata);
outRGB=imageRGB;
for i=1:noX
    for j=1:noY
        if mask(i,j)==0
            outRGB(i,j,:)=[1,1,1];
        end
    end
end