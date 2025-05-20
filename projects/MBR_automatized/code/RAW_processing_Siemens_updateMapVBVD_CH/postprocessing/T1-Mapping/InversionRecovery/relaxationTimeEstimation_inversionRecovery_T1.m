%% Relaxation Time Estimation T1 based on TI
% Ludger Starke
%
% 19-10-2018: Thomas W. Eigentler
%               Adapted to Siemens


clear, close all

%% Preparation
% If not already done, execute T1MappingPreperation and merge Siemens Raw
% data to one Data set and open the data Set here

config.pathname = 'P:\Measurement_7T\20181020_T1T2_IkeaPhantom\';

% 1. select RAW file if not set already in config
if not(isfield(config,'filename')) && isfield(config,'pathname')
    [config.filename, config.pathname, ~] = uigetfile( ...
        {'*.mat','MergedFiles (*.mat)'; ...
         '*.*',  'All Files (*.*)'}, ...
        'Pick a Meas file', config.pathname);
elseif not(isfield(config,'filename'))
    [config.filename, config.pathname, ~] = uigetfile( ...
        {'*.mat','MergedFiles (*.mat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Meas file', 'MergedFile.mat');
end

load([config.pathname config.filename]);

data = 10^8*abs(Temp.combinedImage);
TIs = Temp.TI/1000;

addpath([pwd filesep 'convenienceFunctions'])                                                % Adds the subfolder 'convenienceFunctions' which contains all the used subfunctions


%% Generate signal mask

% close all
% 
% meanImage = mean(data,3);
% 
% backgroundMask = createBackgroundMask(size(meanImage),15);
% measuredMeanNoise = std(meanImage(backgroundMask));
% meanSigma = measuredMeanNoise/0.6551;
% 
% signalMask = abs(meanImage) > 50*meanSigma;
% 
% masksFigure = figure;
% 
%     temp = abs(meanImage);
%     Max = maxN(temp);
%     maskDisplay(:,:,1) = temp/Max/3 + backgroundMask/6;
%     maskDisplay(:,:,2) = temp/Max/3 + signalMask/6;
%     maskDisplay(:,:,3) = temp/Max/3;
% 
%     imshow(maskDisplay)
%     
%     set(masksFigure,'position',[0 0, 900 700]) 
%     set(masksFigure,'PaperPositionMode','Auto') 
% 
%     imwrite(imresize(maskDisplay,3,'nearest'),['SignalMask', '.png'],'png')
    
%% Draw Own Mask

close all

imagesc(abs(data(:,:,1)));
axis image
str=('Choose a region for the T1 calculation')

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
    
%% Get data voxel series

nDataVoxels = sumN(signalMask);

dataVoxels = zeros(nDataVoxels,numel(TIs));

for ii = 1:numel(TIs)
    
    temp = abs(data(:,:,ii));
    dataVoxels(:,ii) = temp(signalMask);
    
end

%% Fit specification

% 2 param fit
func = @(b,x) abs(b(1)*(1 - 2*exp(-b(2)*x)));
opts = optimset('MaxFunEvals',50000, 'MaxIter',10000,'TolFun',10^(-2),'FunValCheck','on');
b0 = [10^(-5), 1/1000];

%% Get T1 map

T1map = zeros(nDataVoxels,1);

for ii = 1:nDataVoxels

    b0 = [max(dataVoxels(ii,:)), 1/1000];
    OLS = @(b) sum((func(b,TIs) - dataVoxels(ii,:)).^2);
    B = fminsearch(OLS, b0, opts);
    
    T1map(ii) = 1/B(2);

end

meanT1 = mean(T1map);

%% Plot fit for example

voxel = round((rand(1)*(nDataVoxels+1)-0.5));

OLS = @(b) sum((func(b,TIs) - dataVoxels(voxel,:)).^2);
B = fminsearch(OLS, b0, opts);
T1 = 1/B(2);

T = linspace(0,max(TIs)*1.1,10000);
fitData = func(B,T);

[~] = standardFigureDefaults();
Fig = figure;

    setupFigure(Fig,[200,200,1000,800])
    load('matlabColors')
    
    plot(TIs,dataVoxels(voxel,:),'+','Color',colors{1},'markersize',12)
    hold on
    plot(T,fitData,'-','Color',colors{2})
    
    xlabel('TI (ms)')
    ylabel('Signal (arbitrary units)')
    
    title(sprintf('Voxel %d, T1 = %.1f ms',voxel, T1)) 

%% Plot T1 map

T1image = zeros(size(meanImage));
T1image(signalMask) = T1map;

figure
imshow(T1image,[])
colormap jet
colorbar

fprintf('\n--- The estimated mean T1 is %.2f ms ---\n\n',meanT1)   