% SIMPLEGRAPPA performs a simple and slow GRAPPA reconstruction
%   GRAPPA reconstruction is taken from the ESPIRiT toolbox 
%
%   USAGE:
%
%   [reskGRAPPA,kSpaceWithRef]=simpleGRAPPA(kSpace,refLines,kernelSize,regularization)
%
%       kSpace = kSpace [Nx,Ny,Nc] without reference lines
%       refLines = [Nx,Nref,Nc], reference lines
%       kernelSize = 2D GRAPPA kernel size

% 2016-10-26, created, Till Huelnhagen

function [reskGRAPPA,kSpaceWithRef]=simpleGRAPPA(kSpace,refLines,kernelSize,regularization)

    addpath(genpath('D:\Matlab\ESPIRiT'));
    
    % calculate acceleration factor
    dataCol=kSpace(floor(end/2),:,1,1,1,1);
    R=round(numel(dataCol)/sum(dataCol==0));
    
    if nargin < 4
        regularization=0.01;
    end

    % insert reference lines into kSpace
    nRefLines=size(refLines,2);
    kSpaceWithRef=kSpace;

    % correct for asymmetric echo if needed
    sizeDiff=size(kSpace,1)-size(refLines,1);
    refLines(end+1:size(kSpace,1),:)=0;
    refLines=circshift(refLines,[sizeDiff 0 0 0 0 0]);

    % match size of refLines and kSpace by zeropadding 
    refLinesFull=refLines;
    refLinesFull(:,end+1:size(kSpace,2),:)=0;
    refLinesFull=circshift(refLinesFull,[0 floor(size(kSpace,2)/2-nRefLines/2)-R+2 0]); % -1 is vital for R=3, otherwise the kSpace is inconsistent which causes artifacts

    % fill missing lines into kSpace
    kSpaceWithRef(refLinesFull~=0)=refLinesFull(refLinesFull~=0);

    % GRAPPA reconstruction
    %disp('Performing GRAPPA reconstruction ')
    reskGRAPPA = GRAPPA(double(kSpaceWithRef),double(refLines),kernelSize,regularization);

