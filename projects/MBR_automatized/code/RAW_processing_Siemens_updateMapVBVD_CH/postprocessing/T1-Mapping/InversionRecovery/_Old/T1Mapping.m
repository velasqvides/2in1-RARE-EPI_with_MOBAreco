function [T1Map,R]=T1Mapping(combinedImage,MrProt,LLcorr)
% This function calculates a T1 map from a 2D Inversion recovery dataset
% by fitting a 3 parameter equation using a Levenberg-Marquardt non linear
% least square fitting routine
%
% INPUT: combinedImage = complex or absolute images with dimentions X,Y,TI
%        MrProt = Struct with MR protocol parameters (i.e. MrProt.TI)
%        LLcorr = perform Look-Locker correction (when using SSFP or FLASH 
%                 readout), 1=yes, 0=no
%
% OUTPUT: T1Map = T1 map
%         R = Pearson correlation coefficient map (optional), describes fit
%             quality
%
% ToDo: The processSiemensRaw function does not extract the TI times
% correctly yet. Once this is fixed, T1 mapping should work.
%
% 2013 Matthias Dieringer: matthias.dieringer@charite.de
% 2018-10-17: Thomas Eigentler: Progress bar deleted

if(nargin<3)
    LLcorr=0;
end
if (LLcorr == 0)
    disp('Look-Locker correction will not be used!')
end

%get inversion times in [s]
t=(MrProt.TI/1E6)';

% take absolute value and scale it to 1
combinedAbsImage=abs(combinedImage);
combinedAbsImage=combinedAbsImage./max(max(max(combinedAbsImage)));

%get data dimensions
sr=size(combinedAbsImage,1);
sp=size(combinedAbsImage,2);
NI=length(t);

%create an erosion mask (thresholding)
[Mask,~,~]=erosion_mask(combinedAbsImage(:,:,1));
for n=1:NI
    combinedAbsImage(:,:,n)=combinedAbsImage(:,:,n).*Mask;
end
close(gcf)

% start T1 Mapping
disp('T1 Mapping started...')
tic
x0 = [0.02,0.5,.1]; % starting variables (1) Offset (2) Mo (3) T1 [s]

%preallocate
T1Map=zeros(sr,sp);
R=zeros(sr,sp);
OS=zeros(sr,sp);
SL=zeros(sr,sp);
tmp1=t;

%create progressbar
% info.title='T1 Mapping progress...';
% info.size=2;
% pb=progbar(info);

%start fitting
for n=1:sr
    for f=1:sp
        if(combinedAbsImage(n,f,1)~=0) %filter noise
            tmp1(:)=combinedAbsImage(n,f,:); %tmp1 contains now all evaluations of 1 pixel
            res = @(x)  abs( x(1) - (x(2)) * (exp(-(t/x(3)))))-tmp1; % Column vector of residuals (prototype fitting function)
            [x,~,~] = LMFnlsq(res,x0); %perform the fitting algorithm
            OS(n,f)=x(1);
            SL(n,f)=x(2);
            T1Map(n,f)=x(3);
            rsq = corrcoef(tmp1(:), abs(x(1) - x(2) * exp(-(t/x(3)))));
            R(n,f)=round(rsq(2)*1000)/1000;
        end
    end
    % progbar(pb,100/sr*n)
end
% close(pb)

% apply Look-Locker correction
if(LLcorr)
    T1Map(n,f)=T1Map(n,f)*(SL(n,f)/OS(n,f)-1);
end
T1Map(T1Map<0 | isnan(T1Map))=0;

disp(['T1 mapping finished! (Duration: ', num2str(floor(toc/60)), ' minutes and ', num2str(round(rem(toc,60))), ' seconds)'])

%display T1-map
figure;
imagesc(round(T1Map*1000))
axis image; axis off
ch=colorbar;
ylabel(ch, 'T_1 [ms]','FontSize',20)
set(gca,'FontSize',12)
title(['Mean T_1 = ', num2str(round(1000*mean2(T1Map(T1Map>0)))),' \pm ', ...
    num2str(round(1000*std2(T1Map(T1Map>0)))),' ms'])

%display fit quality map (optional)
if(nargout==2)
    figure;
    imagesc(R)
    axis image; axis off
    caxis([0 1])
    ch=colorbar;
    ylabel(ch, 'Correlation coefficient','FontSize',20)
    set(gca,'FontSize',12)
    title(['Mean correlation coefficient = ', num2str(mean2(R(R>0))),' \pm ', ...
        num2str(std2(R(R>0)))])
end
end