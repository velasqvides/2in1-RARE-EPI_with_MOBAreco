%% Relaxation Time Estimation T2 based on TE
% Ludger Starke
%
% 19-10-2018: Thomas W. Eigentler
%               Adapted to Siemens


clear, close all

%% Preparation and Data Import

clear config;
config.pathname = 'P:\Measurement_7T\20181020_T1T2_IkeaPhantom';

config.preprocessing.REMOVE_OVERSAMPLING = 1;
processSiemensRaw;

data = 10^8*abs(combinedImage);
TEs = TE*1000;

set(0,'DefaultAxesFontSize', 16)
set(0,'defaultLineMarkerSize', 8)

addpath(genpath([pwd filesep 'convenienceFunctions']))                      %Adds the subfolder 'convenienceFunctions' which contains all the used subfunctions

TEs = TEs(2:end);                                                           %Starting with the second Echo (pre-saturation!!)
data = data(:,:,2:end);

%% Define signal masks
% 
% signalRoiCenter = [94, 87.5];
% signalRoiRadius = 55;
% 
% backgroundRoiSide = 25;
% 
% signalMask = createCircle([size(data,1),size(data,2)],signalRoiCenter,signalRoiRadius);
% backgroundMask = createBackgroundMask([size(data,1),size(data,2)],backgroundRoiSide);
% 
% masksFigure = figure(1);
% 
%     j = 1;     % Which image to show
% 
%     Max = max(max(data(:,:,j)));
%     maskDisplay(:,:,1) = data(:,:,j)/Max/3 + backgroundMask/6;
%     maskDisplay(:,:,2) = data(:,:,j)/Max/3 + signalMask/6;
%     maskDisplay(:,:,3) = data(:,:,j)/Max/3;
% 
%     imshow(maskDisplay)
%     
%     imwrite(imresize(maskDisplay,3,'nearest'),['T2_mask_', num2str(signalRoiRadius,'%d') '.png'],'png')
% 
%     set(masksFigure,'position',[0 0, 900 700]) 
%     set(masksFigure,'PaperPositionMode','Auto') 

%% Draw Own Mask

close all

imagesc(abs(data(:,:,1)));
axis image
str=('Choose a region for the T2 calculation')

meanImage = mean(data,3);
backgroundMask = createBackgroundMask(size(meanImage),15);
measuredMeanNoise = std(meanImage(backgroundMask));
meanSigma = measuredMeanNoise/0.6551;

signalMask = roipoly;

maskFigure = figure;
temp = abs(meanImage);
Max = maxN(temp);
maskDisplay(:,:,1) = temp/Max/3 + backgroundMask/6;
maskDisplay(:,:,2) = temp/Max/3 + signalMask/6;
maskDisplay(:,:,3) = temp/Max/3;

imagesc(maskDisplay)
axis image

%% Get signal and background std estimates    
    
measuredSignal = zeros(size(data,3),1);
correctedSignal = zeros(size(data,3),1);

measuredNoise = zeros(size(data,3),1);

for j=1:length(measuredSignal)
    squeezedData = data(:,:,j);
    measuredSignal(j) = mean(squeezedData(signalMask));
    measuredNoise(j) = std(data(backgroundMask));
end

sigma = mean(measuredNoise)/sqrt(2 - pi/2);

%% Get data voxel series

nDataVoxels = sumN(signalMask);

dataVoxels = zeros(nDataVoxels,numel(TEs));

for ii = 1:numel(TEs)
    
    temp = abs(data(:,:,ii));
    dataVoxels(:,ii) = temp(signalMask);
    
end

%% Bias corrected signal estimates

riceLookupTable = load('riceMeanLookupTable');

for j=1:length(measuredSignal)
    correctedSignal(j) = lookupTableEstimate(measuredSignal(j),sigma,riceLookupTable);
end

signalFigure = figure(2);

    plot(TEs,measuredSignal,'xk','MarkerSize',11,'LineWidth',2)
    hold on
    plot(TEs,correctedSignal,'+b','MarkerSize',11,'LineWidth',2)

    set(signalFigure,'position',[0 0, 800 550]) 
    set(signalFigure,'PaperPositionMode','Auto') 
    saveas(signalFigure,'Signal','jpg') 

    xlim([0, max(TEs)+50])
    
    xlabel('TE (ms)')
    ylabel('Signal (arbitrary units)')
    
    legend('Not Corrected','Corrected','location','NorthEast')


%% Fit

snrThreshold = 0.7;
index = find(correctedSignal < snrThreshold*sigma, 1, 'first');

if isempty(index)
    index = numel(TEs)+1;
end

TEs_new = TEs(1:index-1);
correctedSignal_new = correctedSignal(1:index-1);


func = @(b,x) b(1)*exp(-b(2)*x) + b(3);
opts = optimset('MaxFunEvals',50000, 'MaxIter',10000,'TolFun',10^(-2),'FunValCheck','on');
b0 = [100,1/400,0];
OLS = @(b) sum((func(b,TEs_new) - correctedSignal_new).^2);

T2map = zeros(nDataVoxels,1);

for ii = 1:nDataVoxels

    b0 = [100, 1/400, 0];
    OLS = @(b) sum((func(b,TEs) - dataVoxels(ii,:)).^2);
    B = fminsearch(OLS, b0, opts);
    
    T2map(ii) = 1/B(2);

end

meanT2 = mean(T2map);

% T = 0:max(TEs)*1.2;
% fitData = func(B,T);

%% Plot fit for example

voxel = round((rand(1)*(nDataVoxels+1)-0.5));

OLS = @(b) sum((func(b,TEs) - dataVoxels(voxel,:)).^2);
B = fminsearch(OLS, b0, opts);
TE = 1/B(2);

T = linspace(0,max(TEs)*1.1,10000);
fitData = func(B,T);

[~] = standardFigureDefaults();
Fig = figure;

    setupFigure(Fig,[200,200,1000,800])
    load('matlabColors')
    
    plot(TEs,dataVoxels(voxel,:),'+','Color',colors{1},'markersize',12)
    hold on
    plot(T,fitData,'-','Color',colors{2})
    
    xlabel('TE (ms)')
    ylabel('Signal (arbitrary units)')
    
    title(sprintf('Voxel %d, T1 = %.1f ms',voxel, TE)) 



%% Plot T1 map

T2image = zeros(size(meanImage));
T2image(signalMask) = T2map;

figure
imshow(T2image,[])
colormap jet
colorbar

fprintf('\n--- The estimated T2 is %f ms with %i TEs ---\n\n',meanT2,numel(TEs)+1)    