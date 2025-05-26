clearvars,
% 0. Parameters.
load(fullfile(pwd,'raw_data','info4Reco.mat'));
nSamples = info4Reco.nSamples;
readoutOversampling = info4Reco.readoutOversampling;
nSpokes = info4Reco.nSpokes;

% 1. Load preprocessed data and corrected 2D trajectory
filePath = fullfile(pwd,'processed_data','scanData');
scanData = readcfl(filePath);
filePath = fullfile(pwd,'processed_data','trajectoryCorrected');
trajectoryCorrected = readcfl(filePath);
% to account for oversampling factor of 2
trajectoryCorrected = ...
    bart(sprintf('scale %4.2f',1/readoutOversampling), trajectoryCorrected);

% 2. Get the spokes requiered for reconstruction
nSpokesReco = nSpokes;
assert(nSpokes >= nSpokesReco);
if nSpokes ~= nSpokesReco
    radialKsp = radialKsp(:,:,1:nSpokesReco,:,:);
    trajectoryCorrected = trajectoryCorrected(:,:,1:nSpokesReco,:);
end

% 4. 3D reconstruction
recoMethod = 'gridding'; % 'NUFFT', 'gridding', 'nlinv', 'pics'

bitmask1 = str2double(evalc("bart('bitmask 3')"));
switch recoMethod
    case 'gridding' % density compensation + adjoint NUFFT
        dcf = densityCompRamLak3D(trajectoryCorrected); % squared Ram-Lak filter
        radialKspDC = radialKsp.*dcf;
        clear radialKsp dcf;
        fileName = fullfile(pwd, 'processed_data', 'radialKspDC');
        writecfl(fileName,radialKspDC);
        clear radialKspDC;
        cmd4 = 'bart nufft -a trajCorrScaled radialKspDC imageVolume';
        cmd5 = sprintf('bart rss %i imageVolume imageVolumeRss',bitmask1);
        cmd6 = sprintf('bart resize -c 0 %i 1 %i 2 %i imageVolumeRss imageGridRing',...
            nSamples,nSamples,nSamples);
        cmd7 = sprintf(['rm radialKspDC*cfl radialKspDC*hdr imageVolume*cfl ...' ...
            'imageVolume*hdr imageVolumeRss*cfl imageVolume*hdr']);
        [status,cmdout] = system(strjoin({cmd3, cmd4, cmd5, cmd6, cmd7}, ';'))
        imageGrid = readcfl(fullfile(pwd, 'processed_data', 'imageGridRing'));
    case 'NUFFT' % inverse NUFFT
        fileName = fullfile(pwd, 'processed_data', 'radialKsp');
        writecfl(fileName,radialKsp);
        cmd6 = 'bart nufft -i trajCorrScaled radialKsp imageVolume';
        clear radialKsp;
        cmd7 = sprintf('bart rss %i imageVolume imageVolumeRss',bitmask1);
        cmd8 = sprintf('bart resize -c 0 %i 1 %i 2 %i imageVolumeRss imageNufft',...
            nSamples,nSamples,nSamples);
        cmd9 = sprintf(['rm radialKsp*cfl radialKsp*hdr imageVolume*cfl imageVolume*hdr ...' ...
            'imageVolumeRss*cfl imageVolume*hdr']);
        [status,cmdout] = system(strjoin({cmd1, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9}, ';'))
        imageNufft = readcfl(fullfile(pwd, 'processed_data', 'imageNufft'));
    case 'pics' % parallel imaging+compressed sensing+density compensation as preconditioner
        dcf = calculate2DradialDCF(trajCorrScaled); % normal Ram-Lak filter
        fileName = fullfile(pwd, 'processed_data', 'dcf_file');
        writecfl(fileName,dcf); % this weights have to be passed to pics command in a .cfl file
        clear dcf;
        fileName = fullfile(pwd, 'processed_data', 'radialKsp');
        writecfl(fileName,radialKsp);
        cmd6 = 'bart nufft -i trajCorrScaled radialKsp coilImg';
        mask4 = str2double(evalc("bart('bitmask 0 1 2')"));
        cmd7 = sprintf('bart fft -u %i coilImg cartesianKsp',mask4);
        clear coil_img;
        cmd8 = 'bart ecalib cartesianKsp sensMaps'; % coil sensitivities using ESPIRiT
        mask5 = str2double(evalc("bart('bitmask 0 1 2')"));
        % -R W:%i:0:0.005 l1-wavelet rgularization and lambda=0.005, -i200 iterations,
        % -p dcf_file pattern of weights
        cmd9 = sprintf(['bart pics -S -R W:%i:0:0.013 -i200 -e -p dcf_file -t trajCorrScaled radialKsp sensMaps imageVolume'],mask5);
        clear radialKsp;
        cmd10 = sprintf('bart resize -c 0 %i 1 %i 2 %i imageVolume imagePics',...
            nSamples,nSamples,nSamples);
        [status,cmdout] = system(strjoin({cmd1, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9, cmd10}, ';'));
        imagePics = readcfl(fullfile(pwd, 'processed_data', 'imagePics'));
end
