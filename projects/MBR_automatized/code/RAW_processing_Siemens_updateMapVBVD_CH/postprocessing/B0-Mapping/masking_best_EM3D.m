function [Mask,Perim,Image]=masking_best_EM3D(ref_im)
% This function erodes a mask starting from the edges of the image.
% ToDo: The initial guess should be a little better... Matthias Dieringer

Mask=zeros(size(ref_im));
Perim=zeros(size(ref_im));

Matrix=size(ref_im);
Partitions=size(ref_im,3);

MaxVal=max(max(max(ref_im)));

% thresh=5; MD
thresh=MaxVal/100;
reply=0;
while 1
    Image=ref_im;
    thresh=thresh*(1+reply/10);
    for i=1:Partitions
        Mask(:,:,i)=imerode(imfill(bwareaopen(ref_im(:,:,i)>thresh,8),'holes'),strel('diamond',0));
        Perim(:,:,i)=bwperim(Mask(:,:,i));
        Image(:,:,i)=Image(:,:,i).*(1-Perim(:,:,i))+Perim(:,:,i)*MaxVal;
        Image_plot3D(1,Partitions,i,0,MaxVal,'gray',Image(:,:,i),0)
        
    end
    
    reply = input('To change masking threshold level, enter e.g. +5 or -5, to stop press [Enter] ');
    % prevent window from disappearing
    shg
    reply(isempty(reply))=0;
    if (reply > 0  || reply < 0)
        continue
    else
        break
    end
    
end
end
