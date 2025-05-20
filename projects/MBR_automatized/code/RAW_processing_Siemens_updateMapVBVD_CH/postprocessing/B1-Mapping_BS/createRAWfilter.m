function RAWFilter=createRAWfilter(FilterType,M,N,Strength)
% createRAWfilter(FilterType,M,N,Strength)
%
% Creates a 2D RAW data filter with size MxN, where the smallest
% dimension determines the filter width.
%
% INPUT:
%           FilterType  : 'Gaussian', 'Hanning', 'Hamming', or 'Blackman'
%           MxN         : Filter dimensions along M and N
%           Strength    : Filter slope (optional and for Gaussian only),
%                         can be negative or positive, default value is 3.
% OUTPUT:
%           RAWFilter   : [MxN] sized 2D filter
%
% Example:
% RawFilterGaussian=createRAWfilter('Gaussian',64,64,3);
%
% For informatioon on the performance of some of these filters, look here:
% A.M. Aibinu, et al; MRI Reconstruction Using Discrete Fourier Transform:
% A tutorial; World Academy of Science, Engineering and Technology 18 2008
%
% 2012 Matthias Dieringer
% Matthias.Dieringer@charite.de

if(nargin<4)
    Strength=3;
elseif(nargin==4 && ~strcmp(FilterType,'Gaussian'))
    disp('Warning: "Strength" has no effect with filter types other than "Gaussian"')
end

[X Y]= meshgrid(1:M,1:N);
if(strcmp(FilterType,'Gaussian'))
    %     Sigma=min(N,M)/Strength; %alternatively
    RAWFilter=exp(-(X-M/2).^2./M^2*Strength-(Y-N/2).^2./N^2*Strength);
else
    % create elliptical mask to truncate unwanted side lobes of the filter
    SizeEllipse=min(N,M);
    %ToDo: vectorize
    Mask=ones(N,M);
    for m=1:M
        for n=1:N
            if(sqrt((N/2+n-N)^2+(M/2+m-M)^2)>=SizeEllipse/2)
                Mask(n,m)=0;
            end
        end
    end
    
    if(strcmp(FilterType,'Hanning'))
        RAWFilter=(1-(.5*(1-cos(2*pi*sqrt((X-M/2).^2+(Y-N/2).^2)/(SizeEllipse-1))))).*Mask;
    elseif(strcmp(FilterType,'Hamming'))
        RAWFilter=(1-(0.54-0.46*cos(2*pi*sqrt((X-M/2).^2+(Y-N/2).^2)/(SizeEllipse-1)))).*Mask;
    elseif(strcmp(FilterType,'Blackman'))
        RAWFilter=(1-(0.42-0.5*cos(2*pi*sqrt((X-M/2).^2+(Y-N/2).^2)/(SizeEllipse-1)) + ...
            0.08*cos(4*pi*sqrt((X-M/2).^2+(Y-N/2).^2)/(SizeEllipse-1)))).*Mask;
    else
        disp('Error: Unknown filter type!')
        help createRAWfilter
        return
    end
end