function [Mask,Perim,Image]=erosion_mask(ref_im)
% This function erodes a mask starting from the edges of a 2D image.
% INPUT: ref_im = 2D reference image to apply the mask
%
% OUTPUT: Mask = Binary mask defined by current signal threshold
%         Perim = Binary outline of the mask (perimeter)
%         Image = perimetered image (to simultaneously display image and perim)
%
% simplified 2013 Matthias Dieringer
% matthias.dieringer@charite.de

% preallocate
Mask=zeros(size(ref_im));
Perim=zeros(size(ref_im));

%initial threshold
MaxVal=max(max(ref_im));
thresh=MaxVal/100;
reply=0;
while 1
    Image=ref_im;
    thresh=thresh*(1+reply/10);
    Mask(:,:)=imerode(imfill(bwareaopen(ref_im>thresh,8),'holes'),strel('diamond',0));
    Perim(:,:)=bwperim(Mask);
    Image=Image.*(1-Perim)+Perim*MaxVal;
    imagesc(Image(:,:)); axis image; colormap(gray);
    title('Increase (+5) or decrease (-5) threshold, to stop press [Enter]')
    reply = input('To change masking threshold level, enter e.g. +5 or -5, to stop press [Enter] ');
    % prevent window from disappearing
    shg
    %in case [Enter] is pressed
    reply(isempty(reply))=0;
    if (reply > 0  || reply < 0)
        continue
    else
        break
    end
end
end
