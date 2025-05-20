%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QualityGuidedUnwrap2D implements 2D quality guided path following phase
% unwrapping algorithm.
%
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% This code can easily be extended for 3D phase unwrapping.
% Posted by Bruce Spottiswoode on 22 December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% actual unwrapping function: GuidedFloodFill(...)
% INPUT:
% complex2DImage_in   2dim complex double array containing magnitude and phase information of image
% maskORthreshold     double-> mask is generated with magnitude threshold (unwrapping only if magnitude>threshold)
%                     2dim logical array with same size as complex2DImage
% plotting            0 for no plotting 
%                     1 for output of imagesc-plots
% 
% OUPUT:  
% 2dim real array of unwrapped phases (in rad) 
function output=Unwrap2D_QualityGuided_inputXY(complex2DImage_in, maskORthreshold,plotting,xpoint,ypoint)  
im_mag=abs(complex2DImage_in);
im_phase=angle(complex2DImage_in);
if plotting
    figure; imagesc(im_mag), colormap(gray), axis square, axis off; title('Magnitude-Image'); 
    figure; imagesc(im_phase), colormap(gray), axis square, axis off; title('unwrapped Phase-Image [rad]');
end

if isscalar(maskORthreshold)
    threshold=maskORthreshold;
    
    maxSig=max(im_mag(:));
    minSig=min(im_mag(:));
    im_mag_norm=(im_mag - minSig)/(maxSig-minSig);
  
    im_mask=zeros(size(im_mag));
    im_mask(im_mag_norm>threshold)=1;   
    if plotting
        figure; imagesc(im_mask), colormap(gray), axis square, axis off; title(['magnitude mask, threshold= ',num2str(threshold)]);
    end
elseif size(maskORthreshold)==size(complex2DImage_in)
    
    im_mask=maskORthreshold;
    if plotting
        figure; imagesc(im_mask), colormap(gray), axis square, axis off; title('magnitude mask'); 
    end
else
    
    disp('error in second input-parameter! ')
    return
    
end
%%
im_unwrapped=zeros(size(complex2DImage_in));               %Zero starting matrix for unwrapped phase
adjoin=zeros(size(complex2DImage_in));                     %Zero starting matrix for adjoin matrix
unwrapped_binary=zeros(size(complex2DImage_in));           %Binary image to mark unwrapped pixels

%% Calculate phase quality map
% im_phase_quality=PhaseDerivativeVariance(im_phase,im_mask);   
im_phase_quality=PhaseDerivativeVariance(im_phase);  
%% Identify starting seed point on a phase quality map
% minp=im_phase_quality(2:end-1, 2:end-1); minp=min(minp(:));
% maxp=im_phase_quality(2:end-1, 2:end-1); maxp=max(maxp(:));
% figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off; title('Phase quality map'); 
% uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
% [xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
% disp('cursor position: ');
% disp(round([xpoint,ypoint]));
%% Unwrap
colref=round(xpoint); rowref=round(ypoint);
im_unwrapped(rowref,colref)=im_phase(rowref,colref);                        %Save the unwrapped values
unwrapped_binary(rowref,colref,1)=1;
if im_mask(rowref-1, colref, 1)==1 adjoin(rowref-1, colref, 1)=1; end       %Mark the pixels adjoining the selected point
if im_mask(rowref+1, colref, 1)==1 adjoin(rowref+1, colref, 1)=1; end
if im_mask(rowref, colref-1, 1)==1 adjoin(rowref, colref-1, 1)=1; end
if im_mask(rowref, colref+1, 1)==1 adjoin(rowref, colref+1, 1)=1; end
im_unwrapped=GuidedFloodFill(im_phase, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

if plotting
    figure; imagesc(im_unwrapped), colormap(gray), axis square, axis off; title('Unwrapped phase [rad]'); 
end
output=im_unwrapped;
end
