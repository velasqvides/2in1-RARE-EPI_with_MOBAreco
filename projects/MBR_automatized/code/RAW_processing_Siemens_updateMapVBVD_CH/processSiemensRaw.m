%PROCESSSIEMENSRAW processes Siemens raw files
% Preliminary version with lots of things left ToDo, use results with care
%
% USAGE:
%
%   processSiemensRaw();
%
% Output to workspace:
%
%   config          contains parameters and file information
%   combinedImage   complex reconstructed image (if channel combination is
%   enabled)
%   kSpace          preprocessed kSpace data (e.g. GRAPPPA reconstructed, 
%                   corrected center for asymmetric echo, corrected 
%                   multislice order)
%   MrProt          MR protocol as read from the header
%   newDataDims     contains names of kSpace dimensions
%   sensitivity     sensitivity maps for each channel 
%                   (only for multiecho data that is combined with SVD approach)
%   twix_obj        imacge object to access the raw data directly
%
%   additional scans like refscan, noise etc. (if available in the rawfile)
%   
%   additional data depending on selected postprocessing options
%
%
% PARAMETERS to be set in CONFIG section, or in advance in config struct
% within the calling script
%
% Preprocessing parameters:
%
% config.preprocessing.APPLY_HAMMING_FILTER (default = 0)
%       If true, Hamming filting is applied to kSpace before iFFT to avoid
%       ringing. This will also result in a smoothing of the images.
%
% config.preprocessing.APPLY_NOISE_DECORRELATION (default = 0)
%       If true noise decorrelation of receiver channels is applied for
%       multi channel aquisitions, if a noise scan is found
%
% config.preprocessing.APPLY_TUKEY_FILTER (default = 0)
%       If true, Tukey filting is applied to kSpace before iFFT to avoid
%       ringing. This will also result in a smoothing of the images.
%
% config.preprocessing.COMBINE_CHANNELS (default = 1)
%       If true, channels from multichannel aquisitions are combined. For
%       multi echo data an SVD based approach is used, for single echo
%       data, a magnitude weighted channel combination is used.
%
% config.preprocessing.CORRECT_ASYMMETRIC_ECHO (default = 1)
%       For asymmetric echo acquisitions a part of k-space is not sampled.
%       In order to get a correct reconstruction, the size of k-space has
%       to be adjusted to the actual image dimensions. The unsampled part 
%       of k-space is zerofilled.
%
% config.preprocessing.DISPLAY_NOISE_CORRELATION (default = 1)
%       Displays noise correlation matrix for multi channel acquisitions
%       when decorrelation is selected and noise scan is found
%
% config.preprocessing.DO_AVERAGE (default = 1)
%       Enables or disables averaging for scans acquired with multiple
%       averages.
%
% config.preprocessing.DO_GRAPPA_RECONSTRUCTION (default = 1)
%       Enables or disables reconstruction of unsampled k-space lines for 
%       GRAPPA accelerated acquisitions
%
% config.preprocessing.KEEP_COIL_IMAGES (default = 0)
%       Flag to keep single channel complex images after fft and coil
%       combination
%
% config.preprocessing.ORDER_MULTISLICE (default = 1)     
%       Sometimes the acquisition order is not the same as the geometric 
%       order in multislice scans. The order is saved in the MR protocol. 
%       If is flag is true the slices will be sorted to match the geometric
%       slice order.
%
% config.preprocessing.PHASE_PARTIAL_FOURIER (default = 1)     
%       fills missing part of k-space for phase partial fourier acquisitions
%
% config.preprocessing.REMOVE_OVERSAMPLING (default = 1)
%       Enables or disables the removal of readout oversampling data
%
% config.preprocessing.IGNORE_SEGMENTS (default = 0)
%       Enables or diasables ignoring segment Mdh index 
%
% config.preprocessing.USE_SIMPLE_GRAPPA (default = 0)
%       Enables a simple GRAPPA reconstruction which is working also for
%       higher acceleration factors, but very slow
%
% config.preprocessing.ZEROFILL_PHASE (default = 1)
%       Zerofills k-space for phase resolution < 100% to avoid image 
%       distortion and obtain correct image size in phase encoding 
%       direction for reconstructed images.
% config.preprocessing.IGNORE_SEGMENTS (default = 0);
% config.preprocessing.READ_ADDITIONAL_SCANS (default= 1;

%
% Postprocessing parameters
%
% config.postprocessing.DO_B1_MAPPING (default = 0)
%       Enables calculation of B0, normalized B1+ maps, and absolute B1+
%       maps for various phase based methods such as Bloch-Siegert-Shift,
%       adiabatic Bloch-Siegert Shift, PhiFA-Cup and some more
%       Output: B1Map,B1Map_abs,B0Map
%
% config.postprocessing.DO_T1_MAPPING (default = 0)
%       Enables T1 mapping for scans with multiple inversion times.
%       Output: T1map, R_T1 (Pearson correlation coefficient map, describes
%               fit quality
%       Look Locker correction can be enabled by setting 
%       config.postprocessing.T1Mapping.LLCorr = 1
%
% config.postprocessing.DO_T2_MAPPING (default = 0)
%       Enables T2 mapping for scans with multiple echo times.
%       Output: T2map, R_T2 (Pearson correlation coefficient map, describes
%               fit quality
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGELOG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% 2013-07-10, created by Till Huelnhagen
% 2013-08-23, Till Huelnhagen, added configuration parameters
%             changed reconstruction from ifft to fft
% 2013-08-28, Till Huelnhagen, changed fft back to ifft for mathematical
%             correctness
% 2013-09-05, added first version of 2D SWI postprocessing module
% 2013-10-14, Added T1, T2 mapping to postprocessing (inversion times for
%             T1 mapping still need to be corrected)
% 2013-12-02, Till Huelnhagen, added correction for phase resolution < 100%
% 2014-01-16, new GRAPPA recon implemented (Matthias Dieringer)
% 2014-01-14, last changed by Till Huelnhagen
% 2014-01-14, Till Huelnhagen, implemented check for TXArray data
% 2014-01-15, Matthias Dieringer, added B1 mapping to postprocessing
% 2014-01-29, Olli Kraus, possible to specify a directory with config.pathname
% 2014-03-24, Till Huelnhagen, corrected bug that caused GRAPPA 
%             reconstruction to not work for later echoes for 2D multi-echo 
%             multi-slice and 3D multi-echo scans
% 2014-05-07, Till Huelnhagen, corrected that caused FFT not to be applied for single
%             channel acquisitions
% 2014-09-16, Till Huelnhagen, added T2* mapping support using the fast ARLO approach
% 2015-03-09, Till Huelnhagen, added channel combination support for channel dimensions
%             ColLinChaPhsEco
% 2015-04-29, Till Huelnhagen, added possibility of hamming filtering of
%             k-space to avoid Gibb's ringing
% 2015-07-29, Till Huelnhagen, added flag to keep oversampling
% 2015-10-26, Till Huelnhagen, added option to keep complex single channel
%             coil images when channel combination is used
% 2015-10-27, Till Huelnhagen, added option to apply noise decorrelation of
%             multiple receiver channels, if noise scan is available. Added
%             option to display noise correlation matrix.
% 2015-11-13, Till Huelnhagen, fixed hamming filter calculation for 3D data.
%             Added option to use Tukey filter.
% 2016-08-01, Till Huelnhagen, changed channel combination of phase images 
%             for single echo scans to sum of complex signal, because this
%             provides better phase images than the so far used magnitude
%             based combination. Still this approach is not SNR optimal.
% 2016-10-12, Olli Weinberger, a new B1Mapping function is now used for
%             phase sensitive B1Mapping (Bloch-Siegert, PhiFA-Cup, ...)
% 2016-10-24, Till Huelnhagen, fixed bug that resulted in incorrect image
%             sorting in array after SVD channel combination.
%             Added check if header file does exist already to skip
%             headder extraction in that case.
% 2016-10-26, Till Huelnhagen, added different GRAPPA reconstruction which
%             is also working for higher acceleration factors, but is rather 
%             slow (see config.preprocessing.USE_SIMPLE_GRAPPA)
% 2017-10-26, Till Huelnhagen, added phase partial fourier reconstruction
% 2017-10-27, Till Huelnhagen, changed tukey and hamming filtering to work
%             correctly also for asymmetric k-space
% 2017-07-27, last changed by Till Huelnhagen
% 2019-08-15, Carl Herrmann, added newest version of mapVBVD to support
%             reading data from 3T (VE11) 
%             1) Create function to check for heaxdecima in MrProt and convert
%                them to decimal
%             2) Disabled the writing of header file, since the new
%                twix_map_obj reads the header
%             3) In case of multi-RAID files (VE11) loop over twix_obj cell
%                and write data in arrays
%             4) Disabled the check for TXArray data, because this requires
%                to read header file, this has to be implemented
% 2019-08-19, Carl Herrmann, added reading gradient delay calibration data
%             from twix_obj
% 2019-08-19, Carl Herrmann, added optional parameter for reading raw data 
%             config.preprocessing.IGNORE_SEGMENTS, which enables or
%             disables the ignoring of segment Mdh index, when reading data
%             into arrays. Default is 0.
% 2019-08-19, Carl Herrmann, save gradient delay calibration in struct with
%             arrays for each measured spatial direction
% 2019-08-19, Carl Herrmann, added optional parameter config.preprocessing.
%             READ_ADDITIONAL_SCANS. Enables or disabels reading additional
%             scans from multi-RAID files (VE11), despite the one
%             containing 'image' data. 
% see SVN log for full history

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: These are just default values. External or not deleted previous
%          configurations which are stored in config structure will 
%          override the local settings! This allows the setting of
%          parameters in a calling function.

% initialize config struct if not existing
if not(exist('config','var'))
    config = struct;
    config.preprocessing = struct;
    config.postprocessing = struct;
end
if not(isfield(config, 'preprocessing'))
    config.preprocessing = struct;    
end
if not(isfield(config, 'postprocessing'))
    config.postprocessing = struct;    
end


% check if parameters are set in calling function, if not, set default
% values

% preprocessing parameters
if not(isfield(config.preprocessing,'APPLY_HAMMING_FILTER'))
    config.preprocessing.APPLY_HAMMING_FILTER = 0;
end
if not(isfield(config.preprocessing,'APPLY_NOISE_DECORRELATION'))
    config.preprocessing.APPLY_NOISE_DECORRELATION = 0;
end
if not(isfield(config.preprocessing,'APPLY_TUKEY_FILTER'))
    config.preprocessing.APPLY_TUKEY_FILTER = 0;
end
if not(isfield(config.preprocessing,'COMBINE_CHANNELS'))
    config.preprocessing.COMBINE_CHANNELS = 1;
end
if not(isfield(config.preprocessing,'CORRECT_ASYMMETRIC_ECHO'))
    config.preprocessing.CORRECT_ASYMMETRIC_ECHO = 1;
end
if not(isfield(config.preprocessing,'DISPLAY_NOISE_CORRELATION'))
    config.preprocessing.DISPLAY_NOISE_CORRELATION = 0;
end
if not(isfield(config.preprocessing,'DO_AVERAGE'))
    config.preprocessing.DO_AVERAGE = 1;
end
if not(isfield(config.preprocessing,'DO_GRAPPA_RECONSTRUCTION'))
    config.preprocessing.DO_GRAPPA_RECONSTRUCTION = 1;
end
if not(isfield(config.preprocessing,'KEEP_COIL_IMAGES'))
    config.preprocessing.KEEP_COIL_IMAGES = 0;
end
if not(isfield(config.preprocessing,'PHASE_PARTIAL_FOURIER'))
   config.preprocessing.PHASE_PARTIAL_FOURIER = 1; 
end
if not(isfield(config.preprocessing,'ORDER_MULTISLICE'))
    config.preprocessing.ORDER_MULTISLICE = 1;
end
if not(isfield(config.preprocessing,'USE_SIMPLE_GRAPPA'))
    config.preprocessing.USE_SIMPLE_GRAPPA = 0;
end
if not(isfield(config.preprocessing,'ZEROFILL_PHASE'))
    config.preprocessing.ZEROFILL_PHASE = 1;
end
if not(isfield(config.preprocessing,'REMOVE_OVERSAMPLING'))
    config.preprocessing.REMOVE_OVERSAMPLING = 1;
end
if not(isfield(config.preprocessing,'IGNORE_SEGMENTS'))
    config.preprocessing.IGNORE_SEGMENTS = 0;
end
if not(isfield(config.preprocessing,'READ_ADDITIONAL_SCANS'))
    config.preprocessing.READ_ADDITIONAL_SCANS = 1;
end
if not(isfield(config.preprocessing,'DO_SVDBasedChCombination'))
    config.preprocessing.DO_SVDBasedChCombination = 1;
end
% postprocessing parameters
if not(isfield(config.postprocessing,'DO_T1_MAPPING'))
    config.postprocessing.DO_T1_MAPPING = 0;
end
if not(isfield(config.postprocessing,'DO_T2_MAPPING'))
    config.postprocessing.DO_T2_MAPPING = 0;
end
if not(isfield(config.postprocessing,'DO_T2STAR_MAPPING'))
    config.postprocessing.DO_T2STAR_MAPPING = 0;
end
if not(isfield(config.postprocessing,'DO_B0_MAPPING'))
    config.postprocessing.DO_B0_MAPPING = 0;
end
if not(isfield(config.postprocessing,'B0_MAPPING_savingSTRING'))
    config.postprocessing.B0_MAPPING_savingSTRING=[];
end
if not(isfield(config.postprocessing,'DO_B1_MAPPING'))
    config.postprocessing.DO_B1_MAPPING = 0;
end
if not(isfield(config.postprocessing,'B1_MAPPING_savingSTRING'))
    config.postprocessing.B1_MAPPING_savingSTRING=[];
end
if not(isfield(config.postprocessing,'DO_SWI'))
    config.postprocessing.DO_SWI = 0;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. PARSE RAW DATA AND READ MR PROTOCOL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ToDo: - (correct output order of data in k-space from mapVBVD)
%         - check permutation of additional scans (dimensions might differ)
%         - implement support for TX-array data


% 1. select RAW file if not set already in config
if not(isfield(config,'filename')) && isfield(config,'pathname')
    [config.filename, config.pathname, ~] = uigetfile( ...
        {'*.dat','Meas files VB17 (*.dat)'; ...
        '*.out','Meas files VB17 (*.out)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Meas file', config.pathname);
elseif not(isfield(config,'filename'))
    [config.filename, config.pathname, ~] = uigetfile( ...
        {'*.dat','Meas files VB17 (*.dat)'; ...
        '*.out','Meas files VB17 (*.out)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Meas file', 'meas.dat');
end

if isempty(config.filename), disp('ERROR: No Measurement data are found.'); return; end

filenameRawFile = fullfile(config.pathname, config.filename);

if ~exist(filenameRawFile, 'file') 
    disp(['ERROR: Measurement file does not exist : ' filenameRawFile]); 
    return;
else
    fprintf('\nReading header...\n\n');
end

% This is currently disabled, as mapVBVD reads header in twix_obj, if
% seperate header file is needed, writing it in raw data folder could 
% implemented in future, CH 2019/08/14
% 2. Extract raw header from file if it does not exist yet
% if ~exist(strrep(filenameRawFile,'.dat','.header'),'file')
%     extractHeaderFromSiemensRaw(filenameRawFile);
% else
%     fprintf('\nHeader file already existing. Skipping header extraction.\n\n');
% end

% This is currently disabled, as mapVBVD reads header/protocol data in 
% twix_obj, twix_obj.hdr.MeasYaps has the same structure as MrProt, CH 2019/08/14
% 3. create MR protocol from header
% [~, filen, ~]=fileparts(config.filename);
% [MrProt, ~] = getMrProt(fullfile(config.pathname, sprintf('%s%s',filen,'.header')));
% disp('Saved MR protocol to MrProt');

% disable as this uses header file, implement this check later
% 4. check for TXArray data
% open header file
% fileID=fopen(fullfile(config.pathname, sprintf('%s%s',filen,'.header')));
% inputCells=textscan(fileID,'%s');
% fclose(fileID); clear fileID filen;
% inputCells=inputCells{1,1};
% % search in raw-header for enabled local SAR supervision flag
% VOPS_ENABLED = [];
% VOPS_ENABLED_index = find(~cellfun(@isempty,strfind(inputCells,'<ParamLong."TXASARSupervisionMode">')),1,'first');
% if ~isempty(VOPS_ENABLED_index)
%     VOPS_ENABLED=str2num(inputCells{VOPS_ENABLED_index+2});
% end
% clear inputCells i VOPS_ENABLED_index;
VOPS_ENABLED = [];

% 5. parse image data
fprintf('\nparsing raw data ...\n\n');

% check for software version and 'disable' raw data correction if VE, 
% as it is currently not implemented (see mapVBVD.m)
% lazy software version check (VB or VE?)
fid = fopen(filenameRawFile);
firstInt = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');
status = fclose(fid);
if status ~= 0
    fprintf(['\n ERROR: closing file ' filenameRawFile ' was not successfull \n'])
end
clear status

if and(firstInt < 10000, secondInt <= 64)
    isRawDataCorrect = false;
    disp('Software version: VE (!?)... raw data correction disabled');
else
    isRawDataCorrect = true;
    disp('Software version: VB (!?)... raw data correction enabled');
end
clear firstInt secondInt fid; 

% check for TX_ARRAY data
if isempty(VOPS_ENABLED)
    if config.preprocessing.REMOVE_OVERSAMPLING
        if (config.preprocessing.IGNORE_SEGMENTS)
            if (config.preprocessing.DO_AVERAGE)
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','removeOS','dorawdatacorrect','ignoreSeg');
                else
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','removeOS','ignoreSeg');
                end
            else
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'removeOS','dorawdatacorrect','ignoreSeg');
                else
                    twix_obj=mapVBVD(filenameRawFile,'removeOS','ignoreSeg');
                end
            end
        else
             if (config.preprocessing.DO_AVERAGE)
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','removeOS','dorawdatacorrect');
                else
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','removeOS');
                end
            else
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'removeOS','dorawdatacorrect');
                else
                    twix_obj=mapVBVD(filenameRawFile,'removeOS');
                end
            end
        end
    else
        if (config.preprocessing.IGNORE_SEGMENTS)
            if (config.preprocessing.DO_AVERAGE)
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','dorawdatacorrect','ignoreSeg');
                else
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','ignoreSeg');
                end
            else
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'dorawdatacorrect','ignoreSeg');
                else
                    twix_obj=mapVBVD(filenameRawFile,'ignoreSeg');
                end
            end
        else
            if (config.preprocessing.DO_AVERAGE)
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'doAverage','dorawdatacorrect');
                else
                    twix_obj=mapVBVD(filenameRawFile,'doAverage');
                end
            else
                if isRawDataCorrect
                    twix_obj=mapVBVD(filenameRawFile,'dorawdatacorrect');
                else
                    twix_obj=mapVBVD(filenameRawFile);
                end
            end
        end
    end
else
    fprintf('Detected RAW data with local SAR supervision. This is currently not supported.\n\nRAW processing aborted.\n')
    return;
    %disp('local SAR monitor was enabled! Switching to Siemens tool. WARNING: No additional parameters supported.')
    %twix_obj{config.scanToProcess}.image=readVBVD(filenameRawFile);
end
clear VOPS_ENABLED isRawDataCorrect;

% check if twix_obj is cell array, if not make it a cell array with dim 1
% cell array dim corresponds to number of scans in RAW data file 
% (VD/VE support multiple scans per raw data file)
if ~iscell(twix_obj)
   tmp{1} = twix_obj;
   twix_obj = tmp;
   clear tmp;
end

% make reading the data into arrays compatible with multi-RAID files from
% VE11, by looping over twix_obj
for s=1:size(twix_obj,2)
    % 6. read data into arrays
    fprintf('\nreading data into arrays\n');
    
    if isfield(twix_obj{s},'image')
        kSpace = twix_obj{s}.image();
        IndexImage = s;
        MrProt = twix_obj{s}.hdr.MeasYaps;
        icePara = twix_obj{s}.image.iceParam(7:8,:);
    end
    if (config.preprocessing.READ_ADDITIONAL_SCANS || isfield(twix_obj{s},'image'))
        % check for additional scans and save data into arrays
        if isfield(twix_obj{s},'refscan')
            refscan=twix_obj{s}.refscan();
        end
        if isfield(twix_obj{s},'noise')
            noise=twix_obj{s}.noise();
            IndexNoise = s;
            MrProtNoise = twix_obj{s}.hdr.MeasYaps;
        end
        if isfield(twix_obj{s},'phasecor')
            phasecor=twix_obj{s}.phasecor();
        end
        if isfield(twix_obj{s},'refscanPC')
            refscanPC=twix_obj{s}.refscanPC();
        end
        if isfield(twix_obj{s},'RTfeedback')
            RTfeedback=twix_obj{s}.RTfeedback();
        end
        if isfield(twix_obj{s},'phasestab')
            phasestab=twix_obj{s}.phasestab();
        end
        % read gradient delay calbriation data in array, if data exists 
%         if isfield(twix_obj{s},'gdelcalib_0')
%             gdelcalib_0=twix_obj{s}.gdelcalib_0();
%         end
%         if isfield(twix_obj{s},'gdelcalib_180')
%             gdelcalib_180=twix_obj{s}.gdelcalib_180();
%         end
%         if isfield(twix_obj{s},'gdelcalib_90')
%             gdelcalib_90=twix_obj{s}.gdelcalib_90();
%         end
%         if isfield(twix_obj{s},'gdelcalib_270')
%             gdelcalib_270=twix_obj{s}.gdelcalib_270();
%         end
        % read gradient delay calibration data in struct containing arrays of 
        % data of different spatial directions 
        fn = fieldnames(twix_obj{s});
        fn_gdc = fn( contains(fn,'gdelcalib') == 1 )';
        if ~isempty(fn_gdc)
            for i = 1:size(fn_gdc,2) 
                gdelcalib.(fn_gdc{i}) = twix_obj{s}.(fn_gdc{i})();
            end
        end
        clear i fn*
    end
end
clear s;
% if raw data was recorded using radial trajectory disable channel
% combination and discard coilImages
if MrProt.sKSpace.ucTrajectory==2
    config.preprocessing.COMBINE_CHANNELS = 0;
    config.preprocessing.KEEP_COIL_IMAGES = 0;
end 

% 7. permute data (so that we have [Columns,Lines,Partitions,Slices,Channels,...])
%   ToDo: - check permutation of additional scans (dimensions might differ)
permuteDims=1:numel(twix_obj{IndexImage}.image.dataDims);
if numel(permuteDims > 2)
    permuteDims(1:5)=[1,3,4,5,2];
    kSpace=squeeze(permute(kSpace,permuteDims));
    % permute additional scans
    %     if exist('refscan','var')
    %         refscan=squeeze(permute(refscan,permuteDims));
    %     end
    %     if exist('noise','var')
    %         noise=squeeze(permute(noise,permuteDims));
    %     end
    %     if exist('phasecor','var')
    %         phasecor=squeeze(permute(phasecor,permuteDims));
    %     end
    %     if exist('refscanPC','var')
    %         refscanPC=squeeze(permute(refscanPC,permuteDims));
    %     end
    if exist('gdelcalib_0','var')
        gdelcalib_0=squeeze(permute(gdelcalib_0,permuteDims));
    end
    if exist('gdelcalib_180','var')
        gdelcalib_180=squeeze(permute(gdelcalib_180,permuteDims));
    end
    if exist('gdelcalib_0','var')
        gdelcalib_90=squeeze(permute(gdelcalib_90,permuteDims));
    end
    if exist('gdelcalib_0','var')
        gdelcalib_270=squeeze(permute(gdelcalib_270,permuteDims));
    end

    
    % permute name cell array to match the data
    newDataDims=cell(size(twix_obj{IndexImage}.image.dataDims));
    for i=1:numel(twix_obj{IndexImage}.image.dataDims)
        newDataDims{i}=twix_obj{IndexImage}.image.dataDims{permuteDims(i)};
    end
    
    % squeeze name cell array
    dataSize=twix_obj{IndexImage}.image.dataSize;
    dataSize(:) = dataSize(permuteDims);
    for i=numel(twix_obj{IndexImage}.image.dataDims):-1:1
        % remove singleton dimensions
        if dataSize(i) == 1
            newDataDims(i)=[];
        end
    end
end

% save dimension names in struct
for i=1:numel(newDataDims)
    eval(sprintf('dimNames.%s=%d;',newDataDims{i},i))
end

clear i dataSize filenameRawFile;

%% Noise decorrelation of receiver channels

addpath(genpath('K:\MATLAB\SVN_repository\trunk\RAW_processing\preprocessing')); % Add non linear FFT Code to Matlab search pat

if (config.preprocessing.APPLY_NOISE_DECORRELATION && exist('noise','var') && twix_obj{IndexImage}.image.NCha > 1)
    fprintf('Applying noise decorrelation\n');
    
    % calculate noise correlation matrix
    if config.preprocessing.DISPLAY_NOISE_CORRELATION
        M = size(double(noise'),2);
        Rn = (1/(M-1))*(noise'*noise);
    end
    
    % make sure last dimension in noise array is channels
    chanDim = find(strcmp(twix_obj{IndexNoise}.noise.sqzDims,'Cha'));
    permuteDims=(1:ndims(noise));
    permuteDims(end)=chanDim;
    permuteDims(chanDim)=ndims(noise);
    noise=permute(noise,permuteDims);
    clear IndexNoise
    
    % calculate noise decorrelation matrix from noise array
    dmtx = ismrm_calculate_noise_decorrelation_mtx(noise);
    % apply noise decorrelation matrix to noise array
    noise = ismrm_apply_noise_decorrelation_mtx(noise,dmtx);
    % restore order of array
    %noise=ipermute(noise,permuteDims);
    
    % make sure last dimension in kSpace array is channels
    chanDim = find(strcmp(newDataDims,'Cha'));
    permuteDims=(1:ndims(kSpace));
    permuteDims(end)=chanDim;
    permuteDims(chanDim)=ndims(kSpace);
    kSpace=permute(kSpace,permuteDims);
    % apply noise decorrelation matrix to kSpace and refscan
    kSpace = ismrm_apply_noise_decorrelation_mtx(kSpace,dmtx);
    % restore order of arrays
    kSpace=ipermute(kSpace,permuteDims);
    
    if exist('refscan','var')
        refscan=permute(refscan,permuteDims);
        refscan = ismrm_apply_noise_decorrelation_mtx(refscan,dmtx);
        refscan=ipermute(refscan,permuteDims);
    end
    
    % calculate and display noise correlation matrixafter decorrelation
    if config.preprocessing.DISPLAY_NOISE_CORRELATION
        M = size(double(noise'),2);
        RnDec = (1/(M-1))*(noise'*noise);
        figure;
        subplot(121);imagesc(abs(Rn)); axis image; title('noise correlation matrix');
        subplot(122);imagesc(abs(RnDec)); axis image; colormap(jet); title('noise correlation matrix after decorrelation');
    end
    
    clear dmtx n chanDim permuteDims M Rn RnDec
    
end
    clear permuteDims
    
if (config.preprocessing.APPLY_NOISE_DECORRELATION && ~exist('noise','var'))
    fprintf('\nERROR: Cannot apply noise decorrelation. No noise scan found.\n');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. reconstruct k-space for GRAPPA accelerated data
%    ToDo: - identify flags in MrProt.sPat
%          - R > 2
%          - error in GRAPPA for uneven number of columns
if (hexa2deci(MrProt.sPat.ucPATMode) == 2 && isfield(twix_obj,'refscan') ...
        &&  not(config.preprocessing.DO_GRAPPA_RECONSTRUCTION))
    fprintf('\nWARNING: GRAPPA reconstruction disabled in config\n');
end

%    check for GRAPPA accelerated data
if (hexa2deci(MrProt.sPat.ucPATMode) == 2 && isfield(twix_obj,'refscan') ...
        &&  config.preprocessing.DO_GRAPPA_RECONSTRUCTION)
    
    % start time measurement
    timerGRAPPA = tic;
    
    fprintf('\nDetected 2D GRAPPA acceleration, R = %d\n', MrProt.sPat.lAccelFactPE);
    fprintf('starting GRAPPA reconstruction\n');
    
    % add missing lines to match kSpace size to number of phase encoding lines
    if (size(kSpace,2) < MrProt.sKSpace.lPhaseEncodingLines)
        kSpace(:,end + 1:end + (MrProt.sKSpace.lPhaseEncodingLines - size(kSpace,2)),:) = 0;
    end
    
    % add last 'gap' if needed
    % 2013-08-27, TH, is this necessary? This could make the image size
    %                 larger than it actually is
    if(mod(size(kSpace,2),2)) %avoid odd number of lines
        kSpace(:,end+1,:)=0;
    end
    % do 2D GRAPPA reconstruction in phase encoding direction
    if ndims(kSpace) < 4
        if config.preprocessing.USE_SIMPLE_GRAPPA || MrProt.sPat.lAccelFactPE > 2
            simpleGRAPPA(kSpace,refscan,[3 (MrProt.sPat.lAccelFactPE+1+mod(MrProt.sPat.lAccelFactPE,2))]);
        else
            kSpace=GRAPPA_reko(kSpace,refscan,0);
        end
    else
        % make channel dimension 3rd
        permOrder=1:ndims(kSpace);
        permOrder(3)=dimNames.Cha;
        permOrder(dimNames.Cha)=3;
        kSpace=permute(kSpace,permOrder);
        refscan=permute(refscan,permOrder);
        if config.preprocessing.USE_SIMPLE_GRAPPA || MrProt.sPat.lAccelFactPE > 2
            for i = 1:size(kSpace(:,:,:,:),4)
                fprintf('GRAPPA reconstruction frame %d of %d\n',i,size(kSpace(:,:,:,:),4));
                kSpace(:,:,:,i)=simpleGRAPPA(kSpace(:,:,:,i),refscan(:,:,:,i),[3 (MrProt.sPat.lAccelFactPE+1+mod(MrProt.sPat.lAccelFactPE,2))]);
            end
        else
            for i = 1:size(kSpace(:,:,:,:),4)
                kSpace(:,:,:,i)=GRAPPA_reko(kSpace(:,:,:,i),refscan(:,:,:,i),0);
            end
        end
        kSpace=ipermute(kSpace,permOrder);
        refscan=ipermute(refscan,permOrder);
        clear i permOrder;
    end
    
    fprintf('finished GRAPPA reconstruction\n');
    toc(timerGRAPPA);
    clear timerGRAPPA;
    
elseif (hexa2deci(MrProt.sPat.ucPATMode) ~= 2 || not(isfield(twix_obj,'refscan')))
    fprintf('\nNo GRAPPA acceleration detected\n');
end


% 2. zeropad k-space for asymmetric echo to have the center of k-space in
%    the center of the image and get correct image dimensions
%    ToDo: - Verify cases in which k-space center is offcenter and check
%            suitable preprocessing

if (size(kSpace,1) ~= MrProt.sKSpace.lBaseResolution ...
        && isfield(MrProt.sKSpace, 'AsymmetricEchoAllowed') ...
        && MrProt.sKSpace.ucAsymmetricEchoMode ...
        && not(config.preprocessing.CORRECT_ASYMMETRIC_ECHO))
    fprintf('\nWARNING Asymmetric echo correction disabled in config\n');
end

% check for asymmetric kspace due to asymmetric echo
if (size(kSpace,1) ~= MrProt.sKSpace.lBaseResolution ...
        && isfield(MrProt.sKSpace, 'AsymmetricEchoAllowed') ...
        && MrProt.sKSpace.ucAsymmetricEchoMode ...
        && config.preprocessing.CORRECT_ASYMMETRIC_ECHO)
    
    fprintf('\nDetected asymmetric echo\n');
    fprintf('Adjusting k-space \n');
    
    sizeKSpace = size(kSpace);
    % check if oversampling has already been removed
    if (twix_obj{IndexImage}.image.flagRemoveOS)
        % zerofill kSpace at the top
        kSpaceTmp = zeros([MrProt.sKSpace.lBaseResolution, sizeKSpace(2:end)],'single');
        kSpaceTmp(MrProt.sKSpace.lBaseResolution-sizeKSpace(1)+1:end,:)=kSpace(:,:);
        kSpace = kSpaceTmp;
    else
        % zerofill kSpace at the top
        kSpaceTmp = zeros([MrProt.sKSpace.lBaseResolution * 2, sizeKSpace(2:end)],'single');
        kSpaceTmp(MrProt.sKSpace.lBaseResolution * 2-sizeKSpace(1)+1:end,:)=kSpace(:,:);
        kSpace = kSpaceTmp;
    end
    clear kSpaceTmp sizeKSpace;
end


% 2.1 zeropad k-space for phase resolution < 100% to avoid distortion and
%     get correct image size
%     ToDo: - check phase oversampling

if (MrProt.sKSpace.dPhaseResolution < 1 ...
        && not(config.preprocessing.ZEROFILL_PHASE))
    fprintf('\nWARNING phase zerofilling disabled in config. Reconstructed images will be distorted.\n');
end

centerLin = unique(twix_obj{IndexImage}.image.centerLin);

%% more than one different center lines in encoding
if numel(centerLin) > 1
    fprintf('k-space centerline not consistent, skipping phase zerofilling. Reconstructed images may be distorted.\n');
end

if (MrProt.sKSpace.dPhaseResolution < 1 ...
        && numel(centerLin) == 1 ...
        && config.preprocessing.ZEROFILL_PHASE)
    
    fprintf('\nDetected phase resolution < 100%% \n');
    fprintf('Adjusting k-space \n');
    
    sizeKSpace = size(kSpace);
    % create full size kSpace and copy data of current kSpace
    fullPhaseMatrixSize = MrProt.sKSpace.lBaseResolution*MrProt.sSliceArray.asSlice{1}.dPhaseFOV/MrProt.sSliceArray.asSlice{1}.dReadoutFOV;
    kSpaceTmp = zeros([MrProt.sKSpace.lBaseResolution, fullPhaseMatrixSize ,sizeKSpace(3:end)],'single');
    gap = fullPhaseMatrixSize - size(kSpace,2);
    kSpaceTmp(:, 1:size(kSpace,2),:)=kSpace(:,:,:);
    % calculate new phase center line
    newCenterLin = round(centerLin / size(kSpace,2) * fullPhaseMatrixSize);
    % shift k-space center line to image center
    kSpaceTmp=circshift(kSpaceTmp,double([0,(newCenterLin - centerLin)]));
    kSpace = kSpaceTmp;
    
    clear kSpaceTmp sizeKSpace;
end

clear centerLin;


% 2.2 correct for partial Fourier by filling the not sampled part of kspace
% with the complex conjugate of the opposite positions with regard to the center
% resize kSpace due to partial fourier

if (hexa2deci(MrProt.sKSpace.ucPhasePartialFourier) ~= 16 ...
        && not(config.preprocessing.PHASE_PARTIAL_FOURIER) )
    fprintf('\nWARNING partial Fourier detected, but reconstruction disabled in config. Reconstructed images will be distorted.\n');
end

if (hexa2deci(MrProt.sKSpace.ucPhasePartialFourier) ~= 16 ...
        && config.preprocessing.PHASE_PARTIAL_FOURIER )
    
    fprintf('\nDetected partial Fourier \n');
    fprintf('Adjusting k-space \n');
    %imagesc(log(abs(kSpace(:,:,1)))); axis image; colormap(gray(256)); shg;
    % copy mirrored part of k-space
    fillRange=size(kSpace,2)+1:MrProt.sKSpace.lPhaseEncodingLines;
    % 2D scan
    if twix_obj{IndexImage}.image.NPar < 2
        kSpace(:,fillRange,:,:,:,:,:,:,:,:,:,:)=flipdim(flipdim(conj(kSpace(:,1:numel(fillRange),:,:,:,:,:,:,:,:,:,:)),1),2);
        % 3D scan
    else
        kSpace(:,fillRange,:,:,:,:,:,:,:,:,:,:)=flipdim(flipdim(flipdim(conj(kSpace(:,1:numel(fillRange),:,:,:,:,:,:,:,:,:,:)),1),2),3);
    end
end


% 3. rearrange slices for multislice acquitions to match slice order in DICOM
% ToDo: - verify dimension of kSpace that stores slices (e.g. 3D + multislice)

if config.preprocessing.ORDER_MULTISLICE && twix_obj{IndexImage}.image.NSli > 1
    % check if slices have to be reordered and verify that slices are third
    % dimension of kSpace
    if not(isequal(MrProt.sSliceArray.alSliceAcqOrder, 0:(twix_obj{IndexImage}.image.NSli - 1)))
        fprintf('\nDetected misordered slices:\n');
        fprintf('%d ', cell2mat(MrProt.sSliceArray.alSliceAcqOrder(:)) );
        if strcmp(newDataDims{3},'Sli')
            fprintf('\n\nReordering according to MrProt.sSliceArray.alSliceAcqOrder\n');
            % get order from MrProt
            ord=MrProt.sSliceArray.alSliceAcqOrder;
            % replace empty cell entry with zero, CH 10/21/2019
            if (find(cellfun(@isempty,ord)))
                ord{1,find(cellfun(@isempty,ord))}=0;
                ord = cell2mat(ord);
            else
                ord = cell2mat(ord);
            end
            ord=ord+1;
            % construct array with corrected order
            permuteVec = zeros(size(ord));
            for i=1:twix_obj{IndexImage}.image.NSli
                permuteVec(i) = find(ord==i);
            end
            clear i;
            
            % reorder slices
            kSpace=kSpace(:,:,permuteVec,:,:,:,:,:,:,:,:,:,:,:,:,:);
            fprintf('\nNew slice order:\n');
            fprintf('%d ', cell2mat(MrProt.sSliceArray.alSliceAcqOrder(permuteVec)) );
            fprintf('\n');
        else
            % slice dimension not == 3
            fprintf('\nSlice dimension in kSpace not equal to 3. Slice reordering canceled. Check permutation of kSpace.\n');
        end
    end
    
    % clean up
    clear permuteVec ord;
end

% 4. combine channels, if multiecho, use SVD based approach
%    otherwise use squared magnitude weighted combination
%    combination is done in image space, so fft must be applied first

% ToDo:
% -
if (~config.preprocessing.COMBINE_CHANNELS)
    fprintf('\nSkipping coil channel combination\n');
end
if (config.preprocessing.COMBINE_CHANNELS && twix_obj{IndexImage}.image.NCha > 1 )
    % get position of channels in array
    chanDim = find(strcmp(newDataDims,'Cha'));
    
    if isempty(chanDim)
        disp('ERROR locating channel dimension');
        return;
    end
    
    fprintf('\nStarting channel combination \n');
    
    % apply 2D or 3D fourier transform based on image dimension
    if twix_obj{IndexImage}.image.NPar < 2
        % 2D data, 2D fourier transform along lines and columns
        % Apply hamming window to kSpace to reduce Gibbs ringing
        if config.preprocessing.APPLY_HAMMING_FILTER==1
            fprintf('\nApplying hamming filter to avoid ringing \n');
            % use the actual dimension of the reconstructed image to create
            % the filter and then adjust its size to avoid a displacement
            % between the center of the filter and the center of k-space
            % (e.g. when partial Fourier or asymm. echo are not corrected)
            % which would result in an unequal filtering of positive and
            % negative frequencies
            f = hamming(MrProt.sKSpace.lBaseResolution)*hamming(MrProt.sKSpace.lPhaseEncodingLines)';
            f=f(end-size(kSpace,1)+1:end,1:size(kSpace,2));
            %f = hamming(size(kSpace,1))*hamming(size(kSpace,2))';
            %             % make sure that the filter has the same cutoff frequencies in read and phase direction
            %             sMax=max(MrProt.sKSpace.lBaseResolution,MrProt.sKSpace.lBaseResolution);
            %             f= hamming(sMax)*hamming(sMax)';
            %             if MrProt.sKSpace.lBaseResolution > MrProt.sKSpace.lPhaseEncodingLines
            %                 f=circshift(f,[0 -(MrProt.sKSpace.lBaseResolution-MrProt.sKSpace.lPhaseEncodingLines)/2]);
            %             else
            %                 f=circshift(f,[-(MrProt.sKSpace.lPhaseEncodingLines-MrProt.sKSpace.lBaseResolution)/2, 0]);
            %             end
            %             f=f(end-size(kSpaceFull,1)+1:end,1:size(kSpaceFull,2));
            %             clear sMax
            s = size(kSpace);
            s(1:2)=1;
            kSpace=kSpace.*repmat(f,s);
            clear s f;
        end
        % Apply tukey window to kSpace to reduce Gibbs ringing
        if config.preprocessing.APPLY_TUKEY_FILTER==1
            fprintf('\nApplying tukey filter to avoid ringing \n');
            alpha=0.2;
            %f = tukey(size(kSpace,1),alpha)*tukey(size(kSpace,2),alpha)';
            f = tukey(MrProt.sKSpace.lBaseResolution,alpha)*tukey(MrProt.sKSpace.lPhaseEncodingLines,alpha)';
            f=f(end-size(kSpace,1)+1:end,1:size(kSpace,2));
            s = size(kSpace);
            s(1:2)=1;
            kSpace=kSpace.*repmat(f,s);
            clear s f;
        end
        % apply ifft
        coilImages = ifft(ifftshift(kSpace,1),[],1);
        coilImages = ifft(ifftshift(coilImages,2),[],2);
        coilImages = fftshift(coilImages,1);
        coilImages = fftshift(coilImages,2);
    else
        % 3D data, 3D fourier transform along lines, columns and partitions
        % Apply hamming window to kSpace to reduce Gibbs ringing
        if config.preprocessing.APPLY_HAMMING_FILTER==1
            fprintf('\nApplying hamming filter to avoid ringing \n');
            %f1 = hamming(size(kSpace,1))*hamming(size(kSpace,2))';
            f1 = hamming(MrProt.sKSpace.lBaseResolution)*hamming(MrProt.sKSpace.lPhaseEncodingLines)';
            f1=f1(end-size(kSpace,1)+1:end,1:size(kSpace,2));
            f1 = repmat(f1,[1 1 size(kSpace,3)]);
            f2 = hamming(size(kSpace,3));
            f2 = permute(f2,[2 3 1]);
            f2 = repmat(f2,[size(kSpace,1),size(kSpace,2),1]);
            f = f1.*f2;
            s = size(kSpace);
            s(1:3)=1;
            kSpace=kSpace.*repmat(f,s);
            clear s f1 f2 f;
        end
        % Apply tukey window to kSpace to reduce Gibbs ringing
        if config.preprocessing.APPLY_TUKEY_FILTER==1
            fprintf('\nApplying tukey filter to reduce ringing \n');
            alpha=0.2;
            %f1 = tukey(size(kSpace,1),alpha)*tukey(size(kSpace,2),alpha)';
            f1 = tukey(MrProt.sKSpace.lBaseResolution,alpha)*tukey(MrProt.sKSpace.lPhaseEncodingLines,alpha)';
            f1=f1(end-size(kSpace,1)+1:end,1:size(kSpace,2));
            f1 = repmat(f1,[1 1 size(kSpace,3)]);
            f2 = tukey(size(kSpace,3),alpha);
            f2 = permute(f2,[2 3 1]);
            f2 = repmat(f2,[size(kSpace,1),size(kSpace,2),1]);
            f = f1.*f2;
            s = size(kSpace);
            s(1:3)=1;
            kSpace=kSpace.*repmat(f,s);
            clear s f1 f2 f;
        end
        % apply ifft
        coilImages = ifft(ifftshift(kSpace,1),[],1);
        coilImages = ifft(ifftshift(coilImages,2),[],2);
        coilImages = ifft(ifftshift(coilImages,3),[],3);
        coilImages = fftshift(coilImages,1);
        coilImages = fftshift(coilImages,2);
        coilImages = fftshift(coilImages,3);
    end
    
    % check if multiple echos are available
    if (twix_obj{IndexImage}.image.NEco > 1 && config.preprocessing.DO_SVDBasedChCombination==1)
        fprintf('\nDetected multiecho data, using SVD based channel combination\n');
        
        dimsStr=strcat(newDataDims{:});
        if twix_obj{IndexImage}.image.NPar < 2
            % 2D data,
            % make sure first dimensions are col, row, ch, echo
            if strcmp(dimsStr(1:12),'ColLinChaEco')
                if ndims(coilImages) < 5
                    [sensitivity, combinedImage]=combineSVD2D(coilImages);
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                else
                    siz=size(coilImages);
                    dim = size(coilImages(:,:,:,:,:),5);
                    sensitivity = zeros([siz(1:3),siz(5:end)],'single');
                    combinedImage = zeros([siz(1:2),siz(4:end)],'single');
                    for i = 1:dim;
                        [sensitivity(:,:,:,i), combinedImage(:,:,:,i)]=combineSVD2D(coilImages(:,:,:,:,i));
                    end
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                    % clean up
                    clear dim siz i;
                end
            elseif strcmp(dimsStr(1:15),'ColLinSliChaEco')
                if ndims(coilImages) < 6
                    [sensitivity, combinedImage]=combineSVD(coilImages);
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                else
                    siz=size(coilImages);
                    dim = size(coilImages(:,:,:,:,:,:),6);
                    sensitivity = zeros([siz(1:4),siz(6:end)],'single');
                    combinedImage = zeros([siz(1:3),siz(5:end)],'single');
                    for i = 1:dim;
                        [sensitivity(:,:,:,:,i), combinedImage(:,:,:,:,i)]=combineSVD(coilImages(:,:,:,:,:,i));
                    end
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                    clear dim siz i;
                end
            elseif strcmp(dimsStr(1:15),'ColLinChaPhsEco')
                if ndims(coilImages) < 6
                    [sensitivity, combinedImage]=combineSVD(permute(coilImages,[1 2 4 3 5]));
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                else
                    siz=size(coilImages);
                    dim = size(coilImages(:,:,:,:,:,:),6);
                    sensitivity = zeros([siz(1:4),siz(6:end)],'single');
                    combinedImage = zeros([siz(1:3),siz(5:end)],'single');
                    for i = 1:dim;
                        [sensitivity(:,:,:,:,i), combinedImage(:,:,:,:,i)]=combineSVD(coilImages(:,:,:,:,:,i));
                    end
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                    clear dim siz i;
                end
            elseif strcmp(dimsStr,'ColLinSliChaPhsEco')
                % make sure first dimensions are col, row, sl, ch, echo
                coilImages=permute(coilImages,[1 2 3 4 6 5]);
                siz=size(coilImages);
                dim = size(coilImages(:,:,:,:,:,:),6);
                sensitivity = zeros([siz(1:4),siz(6:end)],'single');
                combinedImage = zeros([siz(1:3),siz(5:end)],'single');
                for i = 1:dim;
                    [sensitivity(:,:,:,:,i), combinedImage(:,:,:,:,i)]=combineSVD(coilImages(:,:,:,:,:,i));
                end
                % make echo dimensions last in array
                combinedImage=permute(single(combinedImage),[1 2 3 5 4]);
                coilImages=permute(coilImages,[1 2 3 4 6 5]);
                clear dim siz i;
            else
                fprintf('\nUnknown order of kSpace dimensions. Channel combination skipped.\n');
                % data is ordered differently
                % ToDo!
            end
        else
            % 3D data
            % make sure first dimensions are col, row, par, ch, echo
            if strcmp(dimsStr(1:15),'ColLinParChaEco')
                if ndims(coilImages) < 6
                    [sensitivity, combinedImage]=combineSVD(coilImages);
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                else
                    siz=size(coilImages);
                    dim = size(coilImages(:,:,:,:,:,:),6);
                    sensitivity = zeros([siz(1:4),siz(6:end)],'single');
                    combinedImage = zeros([siz(1:3),siz(5:end)],'single');
                    for i = 1:dim;
                        [sensitivity(:,:,:,:,i), combinedImage(:,:,:,:,i)]=combineSVD(coilImages(:,:,:,:,:,dim));
                        % use single precision to save memory
                    end
                    sensitivity=single(sensitivity);
                    combinedImage=single(combinedImage);
                    clear dim siz i;
                end
            else
                fprintf('\nUnknown order of kSpace dimensions. Channel combination skipped.\n');
                % data is ordered differently
                % ToDo!
            end
        end
        
    else
        % no multiecho acquisition: Use SoS for magnitude and complex sum
        % for phase. Complex sum does not provide optimal SNR, but better
        % phase than magnitude weighted combination
        fprintf('\nNo multiecho, using complex sum of single channel images for phase and SoS for magnitude\n');
        phaseCombined = angle(squeeze(sum(coilImages.*1,chanDim)));
        magnCombined=single(sqrt(squeeze(sum(abs(coilImages).^2,chanDim)/size(coilImages,chanDim))));
                
        % combine magnitude and phase in complex image
        combinedImage = magnCombined.*exp(1i*phaseCombined);
    end
    
    clear chanDim dimsStr magnCombined phaseCombined;
    
    fprintf('\nFinished coil combination\n');
    
    % apply fft also for single channel datasets
else
    % apply 2D or 3D fourier transform based on image dimension
    if twix_obj{IndexImage}.image.NPar < 2
        % 2D data, 2D fourier transform along lines and columns
        coilImages = ifft(ifftshift(kSpace,1),[],1);
        coilImages = ifft(ifftshift(coilImages,2),[],2);
        coilImages = fftshift(coilImages,1);
        coilImages = fftshift(coilImages,2);
    else
        % 3D data, 3D fourier transform along lines, columns and partitions
        coilImages = ifft(ifftshift(kSpace,1),[],1);
        coilImages = ifft(ifftshift(coilImages,2),[],2);
        coilImages = ifft(ifftshift(coilImages,3),[],3);
        coilImages = fftshift(coilImages,1);
        coilImages = fftshift(coilImages,2);
        coilImages = fftshift(coilImages,3);
    end
    % use variable name combinedImage in single channel case, otherwise keep coilImages
    if twix_obj{IndexImage}.image.NCha < 2
        combinedImage = coilImages;
    end
end
if config.preprocessing.KEEP_COIL_IMAGES == 0
    clear coilImages
end
fprintf('\nRAW preprocessing finished\n');  


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG

% imagesc(angle(combinedImage(:,:,9))); axis image; axis off; colormap gray;

% imagesc((magn(:,:,1))); axis image; axis off; colormap gray;

%%%%
% figure;
% subplot(2,2,1)
% imagesc((log(abs(kSpace(:,:,3,1))))); axis image;
% title('log(asym k-space)')
% 
% subplot(2,2,2)
% imagesc(abs(ifftshift(ifft2(ifftshift(kSpace(:,:,3,1)))))); axis image;
% title('magnitude image')
% 
% subplot(2,2,3)
% imagesc((log(abs(kSpaceTest(:,:,3,1))))); axis image;
% title('log(zeropadded k-space)')
% 
% subplot(2,2,4)
% imagesc(abs(ifftshift(ifft2(ifftshift(kSpaceTest(:,:,3,1)))))); axis image;
% title('magnitude image')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. POSTPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the desired postprocessing you like
% ToDo: - Add modules and checks to make sure that the data is suitable for
%         the respective module (e.g. multiecho for B0)


% Provide TE in ms
%TE = MrProt.alTE{1:MrProt.lContrasts}/1e6;
TE = cell2mat(MrProt.alTE);
TE = TE(1:MrProt.lContrasts)/1e6;

% ----------------------------------------------------------------------- %
%   T1 mapping
%       calculates a T1 map. For details see script in subfolder 
%       T1T2Mapping
%   ToDo: - check TI times in MrProt
%         - verify correct input data
% ----------------------------------------------------------------------- %
if (config.postprocessing.DO_T1_MAPPING)
    fprintf('\nT1 mapping\n');  
    if not(isfield(config.postprocessing,'T1Mapping'))
        config.postprocessing.T1Mapping = struct;
    end
    if not(isfield(config.postprocessing.T1Mapping,'LLcorr'))
        config.postprocessing.T1Mapping.LLcorr = 0;
    end
    [T1Map,R_T1]=T1Mapping(combinedImage,MrProt,LLcorr);
end


% ----------------------------------------------------------------------- %
%   T2 mapping
%       calculates a T2 or T2* map respectively. For details see script in 
%       subfolder T1T2Mapping
%   ToDo: - verify correct input data
% ----------------------------------------------------------------------- %
if (config.postprocessing.DO_T2_MAPPING)
    fprintf('\nT2 mapping\n'); 
    [T2Map,R_T2]=T2Mapping(combinedImage,MrProt);
end

% ----------------------------------------------------------------------- %
%   T2* mapping
%       calculates T2* map using ARLO. For details see script in 
%       subfolder T2starMapping
%   ToDo: - verify correct input data
%         - verify equidistant echoes
% ----------------------------------------------------------------------- %
if (config.postprocessing.DO_T2STAR_MAPPING)
    fprintf('\nPerforming log linear T2* mapping with noise cutoff. Result saved in t2starFit struct\n');
    %[T2starMap]=T2StarMappingARLO(combinedImage,MrProt,newDataDims);
    % set noise level to lower 5% of magnitude unless specified externally
    if not(isfield(config.postprocessing,'t2StarMapping'))        
        config.postprocessing.t2StarMapping.noiseLevel = 0.02*max(abs(combinedImage(:)));
    elseif not(isfield(config.postprocessing.t2StarMapping,'noiseLevel'))
        config.postprocessing.t2StarMapping.noiseLevel = 0.02*max(abs(combinedImage(:)));
    end    
    % check number of echoes
    if MrProt.lContrasts < 2
        fprintf('\WARNING: Insufficient number of echoes for T2* fitting. This requires at least 2 echoes\n');
    else
        [t2starFit.t2starMap t2starFit.m0 t2starFit.rsquare t2starFit.adjRsquare t2starFit.rmse t2starFit.sse t2starFit.cropNum t2starFit.t2starStd] = logLinLST2starFitNoiseCut(double(abs(combinedImage)), TE, config.postprocessing.t2StarMapping.noiseLevel);
        fprintf('\nFinished T2* fitting.\n');
    end
end

% ----------------------------------------------------------------------- %
%   B0 mapping
%       This function calculates B0 based on multi-echo GRE sequences
%   ToDo: - check 3D datasets (processSiemensRaw did not reconstruct images
%   properly)
% ----------------------------------------------------------------------- %
if(config.postprocessing.DO_B0_MAPPING)
    fprintf('\nB0+ mapping\n'); 
    [B0Map_Hz,detailedData]=B0Mapping_20140923(combinedImage,MrProt,newDataDims,config);
end

% ----------------------------------------------------------------------- %
%   B1 mapping
%       This function calculates B0, normalized B1+ maps, and absolute B1+
%       maps for various phase based methods such as Bloch-Siegert-Shift,
%       adiabatic Bloch-Siegert Shift, PhiFA-Cup and some more
%   ToDo: - check 3D datasets (processSiemensRaw did not reconstruct images
%   properly)
% ----------------------------------------------------------------------- %
if(config.postprocessing.DO_B1_MAPPING)
    fprintf('\nB1+ mapping\n'); 
    
    %[B1Map_rel,B1Map_uT,B0Map_Hz,detailedData] = myB1Mapping_20141006(combinedImage,MrProt,newDataDims,config);
    %[B1Map_rel,B1Map_uT,B0Map_Hz,detailedData] = B1Mapping_20150521(combinedImage,MrProt,newDataDims,config);
    [Maps,detailedB1postprocessingData] = B1Mapping_20161012(combinedImage,MrProt,newDataDims,config);
    
end

% ----------------------------------------------------------------------- %
%   SWI processing
%       creates a magnitude image, that is enhanced by a phase image
%       according to Haacke et al. (2004) Magn Reson Med 52(3):612.
%       no minimum intensity projection is carried out!
% ----------------------------------------------------------------------- %
%   ToDo: - add MIP as option
%         - verify image dimensions column, line, (echo)
%         - multislice and 3D data
%         - cases with more dimensions (e.g. CINE)

% if (config.postprocessing.DO_SWI && not(exist('combinedImage','var')))
%     fprintf('\nERROR: Unable to do SWI processing. No combined image data found. Check multichannel combination\n');
% end
% 
% if (config.postprocessing.DO_SWI && exist('combinedImage','var'))
%     fprintf('\nstarting SWI processing\n');
%     
%     % check if parameters are already set. If not, set default values
%     if not(isfield(config.postprocessing,'SWI'))
%         config.postprocessing.SWI = struct;
%     end
%     if not(isfield(config.postprocessing.SWI,'maskExponent'))
%         config.postprocessing.SWI.maskExponent = 4;
%     end
%     if not(isfield(config.postprocessing.SWI,'phaseMaskType'))
%         config.postprocessing.SWI.phaseMaskType = 'negative';
%     end
%     if not(isfield(config.postprocessing.SWI,'filterStrength'))
%         config.postprocessing.SWI.filterStrength = 0.05;
%     end
%     if not(isfield(config.postprocessing.SWI,'phaseEcho'))
%         config.postprocessing.SWI.phaseEcho = twix_obj{IndexImage}.image.NEco;
%     end
%     if not(isfield(config.postprocessing.SWI,'magnitudeEcho'))
%         config.postprocessing.SWI.magnitudeEcho = 1;
%     end
% %     
% %     % check for 2D or 3D data
% %     if twix_obj{IndexImage}.image.NPar < 2
% %         % 2D data
% %         % highpassfilter data
% %         hpFilteredImage = homodyne2D(squeeze(combinedImage(:,:,config.postprocessing.SWI.phaseEcho,1,1,1,1,1,1,1,1,1,1,1,1)), config.postprocessing.SWI.filterStrength);
% %         
% %         % create phase mask
% %         hpPhase = angle(hpFilteredImage);
% %         phaseMask = createSWIPhaseMask(hpPhase, config.postprocessing.SWI.phaseMaskType);
% %     
% %         % apply phase mask
% %         SWI = phaseMask.^config.postprocessing.SWI.maskExponent.*squeeze(abs(combinedImage(:,:,config.postprocessing.SWI.magnitudeEcho,1,1,1,1,1,1,1,1,1,1,1,1)));
% %         
% %     else
% %         % 3D data
% %         
% %     end
% %     
% %     % clean up
% %     clear phaseMask hpPhase hpFilteredImage;
% %     
% %     fprintf('\nSWI processing finished\n');
% end
clear IndexImage dimNames