function [dDcm, dDcmInfo,pathname, filenames] = readDcmFolder(varargin)
%*************************************************************************
% imports Dicom images and the corresponding header of one or multiple DCM
% files. The filename can be entered as a varargin or if nothing is
% entered, then a UI will ask for the files. Multiselection is possible.
%
%
%   Modificatons: 
%   20170117: changed dicominfo struct to cell
%
% INPUT:
%   path(s) of DCM files (optional). If nothing is entered, then a UI will
%   ask for the files. Multiselect is possible
%
% OUTPUT:
%   dDcm:       Dicoms
%   dDcmInfo:   headers of the Dicoms
%
% OPTIONS:      DIRNAME     if given, then function does not ask for folder
%
% EXAMPLE:
%---------
%*************************************************************************

propdef = struct;
propdef.VERBOSE = true;
propdef.SORT    = true;
propdef.FILEINCREMENT = 1;
propdef.FILESTARTNO = 1;

prop    = catstruct(propdef,parseVariableInputs(varargin));

%if nargin==0
%if no path is defined as function argument.
if(isfield(prop,'DIRNAME'))
    pathname = prop.DIRNAME;
else
    pathname = uigetdir();
end

dcmFilesDcm = dir(fullfile(pathname,'*.dcm'));
dcmFilesDcm2 = dir(fullfile(pathname,'*.DCM'));
dcmFilesIma = dir(fullfile(pathname,'*.ima'));
dcmFilesIma2 = dir(fullfile(pathname,'*.IMA'));
dcmFiles = catstruct(dcmFilesDcm,dcmFilesIma);
dcmFiles = catstruct(dcmFiles,dcmFilesDcm2);
dcmFiles = catstruct(dcmFiles,dcmFilesIma2);
%this code works only for DCM OR IMA (if both are in the same folder, it
%will fail!)  

lNoOfFiles = size(dcmFiles,1); 

%Anzahl der files ermitteln
%filename = cellstr(filename);
%lNoOfFiles = size(filename,2); 
%filename = circshift(filename,[0 0]);

olddir = cd;
cd(pathname);

%Dicoms einlesen:
if(propdef.VERBOSE)
    fprintf('\n\nReading File #\n\n\n\n\n\n\n\n');
end
lI = 1;
for lL=prop.FILESTARTNO:prop.FILEINCREMENT:lNoOfFiles
   dDcmTmp(:,:,lI) = dicomread(char(dcmFiles(lL).name));
   dDcmInfoTmp{lI} = dicominfo(char(dcmFiles(lL).name)); 
   dSliceLocationTmp(lI) = dDcmInfoTmp{lI}.SliceLocation;
   if(propdef.VERBOSE)
    fprintf('\b\b\b\b\b\b\b\b\b%04d/%04d',lL,lNoOfFiles);
   end
   lI = lI+1;
end
if(propdef.VERBOSE)
    fprintf('\nDone!\n');
end
[sortvalue,sortindex] = sort(dSliceLocationTmp);
%for lL=sortindex
if(prop.SORT)
    dDcm = dDcmTmp(:,:,sortindex);
    dDcmInfo = dDcmInfoTmp{sortindex};  
    dSliceLocation = dSliceLocationTmp(sortindex);
    filenames = dcmFiles(sortindex);
else
    dDcm = dDcmTmp(:,:,:);
    dDcmInfo = dDcmInfoTmp;  
    dSliceLocation = dSliceLocationTmp(:);
    filenames = dcmFiles;
    
end
cd(olddir);

% else
%     
%     lNoOfFiles = nargin;
%     
%     for lL=1:lNoOfFiles
%        dDcm(:,:,lL) = dicomread(char(varargin(lL)));
%        dDcmInfo(lL) = dicominfo(char(varargin(lL))); 
%     end
% 
% end

end