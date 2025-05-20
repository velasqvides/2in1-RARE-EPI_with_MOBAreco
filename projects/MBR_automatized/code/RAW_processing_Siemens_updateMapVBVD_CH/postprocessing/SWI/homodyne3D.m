function [ HPFIm, LPFIm, sig ] = homodyne3D( image, strength )
%HOMODYNE3D applies a homodyne highpassfilter to a complex image input by
%   complex dividing the image by a lowpassfiltered version of it. 4D data
%   are supported as inputs, while the filtering will only be done in the
%   first and second dimension
%
%   USAGE:
%
%       [ HPFIm, LPFIm, sig ] = homodyne3D( image, strength )
%
%   INPUT:
%
%       image       complex valued input image
%       strength    filter strength (the 3dB cutoff (pixels measured from 
%                   the center) as a percent of Kx, e.g. 0.02)
%
%   OUTPUT:
%
%       HPFIm       Highpassfiltered complex image
%       LPFIm       Lowpassfiltered complex image
%       sig         standard deviation of the lowpass filter
%
%
%   05-09-2013, created by Till Huelnhagen, adapted from Grant Kaijuin Yang
%               (Bachelor Thesis, Ohio State University)
%   05-09-2013, last changed by Till Huelnhagen

% Strength gives the strength of the HPF
dimension=size(image);
%check for data sets with only one echo
if size(dimension)<4
    dimension=[dimension,1];
end

%set filter strength (the 3dB cutoff point (pixels measured from the center) as a percent of Kx or Ky)
x1=dimension(1)/2*strength;
sig=sqrt(-.5*x1^2/log10(1/sqrt(2)));
sigky=sig;
sigkx=sig;
%create filter
h=[dimension(1),dimension(2)];

% ToDo: replace by fspecial to avoid dimension/2 which causes problems for
%       odd number of rows or columns
for x=1:dimension(1)
    parfor y=1:dimension(2)
        h(x,y)=exp(-.5*((x-dimension(1)/2)^2/sigkx^2+(y-dimension(2)/2)^2/sigky^2));
    end
end
h=circshift(h,[1 1]); % shift center of h to center of k-space, TH, 2016-08-10

HPFIm=zeros(dimension);
LPFIm=zeros(dimension);
k=0;
for echo=1:dimension(4)
    parfor z=1:dimension(3)
        tmp=fft2(image(:,:,z,echo));
        tmp=fftshift(tmp);        
        tmp=tmp.*h;
        tmp=fftshift(tmp);
        LPFIm(:,:,z,echo)=ifft2(tmp);
        HPFIm(:,:,z,echo)=image(:,:,z,echo)./LPFIm(:,:,z,echo);
    end
    k=k+1;
    fprintf(['Finished Echo ' num2str(k) '\n']);
end

clear tmp;
fprintf('finished homodyne filtering \n');

end