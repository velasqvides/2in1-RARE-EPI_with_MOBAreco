% This example opens a VB 17 RAW file with an Inversion recovery scan
% (GRAPPA R=2, 24ch head coil) calculates the full k-space using GRAPPA
% algorithm, and finally reconstructs the magnitude image
%
% 06-2013 Matthias Dieringer

%start from clear workspace
clear all; close all; clc

% pick and chose the RAW file provided in this folder
[filename, pathname, filterindex] = uigetfile( ...
    {'*.dat','Meas files VB17 (*.dat)'; ...
    '*.out','Meas files VB17 (*.out)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a Meas file', 'meas.dat');

if isempty(filename), disp('No Measurement data are found.'); return; end
dat_fn = fullfile(pathname, filename);

if ~exist(dat_fn, 'file') 
    disp(['Measurement file does not exist : ' dat_fn]); 
    return;
else
    disp('Reading file...');
end

% read RAW image data and reference scan
twix_obj=mapVBVD(dat_fn,'doAverage','removeOS','doRawDataCorrect');
RAWdata = twix_obj.image{''};
refscan=twix_obj.refscan{''};

% permute data (so that we have [Freq,Phase,Coils])
RAWdata=permute(RAWdata,[1,3,2]);
refscan=permute(refscan,[1,3,2]);

% add last 'gap' if needed
if(mod(size(RAWdata,2),2)) %avoid odd number of lines
    RAWdata(:,end+1,:)=0;
end

% look at the undersampled data
figure
imagesc(log(abs(RAWdata(:,:,1))))
title('Undersampled RAW data')
xlabel('Phase encoding')
ylabel('Freq encoding')

% perform GRAPPA reconstruction
tic
full_kspace=Grappa2D(RAWdata,refscan);
% (you can pass over reduction factor manually either)
% full_kspace=Grappa2D(RAWdata,refscan,2);
toc

% look at the GRAPPA data of coil 1
figure
imagesc(log(abs(full_kspace(:,:,1))))
title('GRAPPA reconstructed k-space')
xlabel('Phase encoding')
ylabel('Freq encoding')

% reconstruct images, apply sum of squares, and extract magnitude
image=sqrt(sum(abs(ifftshift(ifft2(ifftshift(full_kspace)))).^2,3)/size(RAWdata,3));

% display final image
figure
imagesc(image)
axis image
axis off
colormap gray