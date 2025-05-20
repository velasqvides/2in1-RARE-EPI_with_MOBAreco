function [T2Map,R]=T2Mapping(combinedImage,MrProt)
% This function calculates a T2 map from 2D Spin Echo Multi Contrast (SEMC)
% MR data by a Levenberg-Marquardt fitting routine using a 3 parameter fit
%
% INPUT: combinedImage = complex or absolute images with dimentions X,Y,TE
%        MrProt = Struct with MR protocol parameters (i.e. MrProt.TE)
%
% OUTPUT: T2Map = T2 map
%         R = Pearson correlation coefficient map (optional), describes fit
%             quality
%
% 2013 Matthias Dieringer
% matthias.dieringer@charite.de
% 2013-10-15, modified by Till Huelnhagen
% 2018-10-17: Thomas Eigentler: Progress bar deleted

%get echo times in [s]
t=(MrProt.TE/1E6)';

% use only as many TE's as there are echos in the acquisition, TH
t=t(1:MrProt.Contrasts);

% discard 1st data point, since it only contains the primary echo
t=circshift(t,[-1,0]); % changed dimensions, TH
t(end)=[];
combinedAbsImage=abs(circshift(combinedImage,[0,0,-1]));
combinedAbsImage(:,:,end)=[];

%get data dimensions
sr=size(combinedImage,1);
sp=size(combinedImage,2);
NI=length(t);

% scale data up to 1 (for more standardized data input)
combinedAbsImage=combinedAbsImage./max(max(max(combinedAbsImage)));

%create a mask
[Mask,~,~]=erosion_mask(combinedAbsImage(:,:,1));
for n=1:NI
    combinedAbsImage(:,:,n)=combinedAbsImage(:,:,n).*Mask;
end
close(gcf)

% start T2 Mapping
disp('T2 Mapping started...')
tic
x0 = [0.02,0.5,.1]; % starting variables (1) Offset (2) Mo (3) T2 [s]

%preallocate
T2Map=zeros(sr,sp);
R=zeros(sr,sp);
tmp1=t;

%create progressbar ---> Deleted TE (Not possible to start Toolbar)
% info.title='T2 Mapping progress...';
% pb=progbar(info);

%start fitting
for n=1:sr
    for f=1:sp
        if(combinedAbsImage(n,f,1)~=0) %filter noise
            tmp1(:)=combinedAbsImage(n,f,:); %tmp1 contains now all evaluations of 1 pixel
            res = @(x)  x(1) + x(2) * (exp(-(t/x(3))))-tmp1; % Column vector of residuals (prototype fitting function)
            [x,~,~] = LMFnlsq(res,x0); %perform the fitting algorithm
            T2Map(n,f)=x(3);
            rsq = corrcoef(tmp1(:),x(1) + x(2) * (exp(-(t/x(3)))));
            R(n,f)=round(rsq(2)*1000)/1000;
        end
    end
    % progbar(pb,100/sr*n)
end
% close(pb)

T2Map(T2Map<0 | isnan(T2Map))=0;

disp(['T2 mapping finished! (Duration: ', num2str(floor(toc/60)), ' minutes and ', num2str(round(rem(toc,60))), ' seconds)'])

%display T2-map
figure;
imagesc(round(T2Map*1000))
axis image; axis off
ch=colorbar;
ylabel(ch, 'T_2 [ms]','FontSize',20)
set(gca,'FontSize',12)
title(['Mean T_2 = ', num2str(round(1000*mean2(T2Map(T2Map>0)))),' \pm ', ...
    num2str(round(1000*std2(T2Map(T2Map>0)))),' ms'])
 
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