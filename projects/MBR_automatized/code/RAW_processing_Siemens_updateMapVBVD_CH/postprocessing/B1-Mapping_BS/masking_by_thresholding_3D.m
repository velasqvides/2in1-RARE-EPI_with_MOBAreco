function [Mask,Image,thresh]=masking_by_thresholding_3D(mag_in)
% This function erodes a mask starting from the edges of the image.
% ToDo: The initial guess should be a little better... Matthias Dieringer

Mask=zeros(size(mag_in));

%Matrix=size(mag_in);
Partitions=size(mag_in,3);

MaxVal=max(mag_in(:));

% thresh=5; MD
thresh=MaxVal/100;
reply=0;
while 1
    Image=mag_in;
    thresh=thresh*(1+reply/10);
    for i=1:Partitions
        Mask(:,:,i)=logical(mag_in(:,:,i)>thresh);
        Image(:,:,i)=mag_in(:,:,i).*Mask(:,:,i);
        Image_plot3D(1,Partitions,i,0,MaxVal,'jet',Image(:,:,i),0,1)        
    end
    %fprintf('Current threshold: %2.2f %%\n',100*thresh/MaxVal);
    reply = input('increase/decrease threshold with a +/- number, to stop press [Enter] ');
    shg % prevent window from disappearing
    reply(isempty(reply))=0;
    if (reply > 0  || reply < 0)
        continue
    else
        break
    end
    
end
end