function [initMaps] = createInitMpasForT2star(kSpace, traj, TEs, sens, R2,lastSynthetic,  whichM0)
bitmask2 = str2double(evalc("bart('bitmask 0 1')"));
% B0_ksp = bart(sprintf('fft -u %i', bitmask2), initB0_fessler);
nSamples = size(kSpace, 2);
nSlices = size(kSpace,7);
initMaps = zeros(nSamples,nSamples,1,1,1,1,3,nSlices);
% at the moment T2 MOBA works for bartv07 only
functionDir = fileparts(mfilename('fullpath'));
bartPath = fullfile(functionDir,'../../../../tools/bart_v07/bart');
run(fullfile(bartPath, 'startup.m'));
for slice=1:nSlices
    TE_e = TEs(:,:,:,:,:,1:3);
    traj_1 = traj(:,:,:,:,:,1:3,slice);
    kdat_1 = kSpace(:,:,:,:,:,1:3,slice);
    R_m0_1fe = bart(sprintf('moba -g -O -G -m0 -i4 -o%f --fat_spec_0 -R3 -t',1),traj_1, kdat_1, TE_e);
    initB0 = bart('extract 6 2 3', R_m0_1fe);

    switch whichM0
        case 0
            initW = ones(nSamples,nSamples,1);
            % initB0 = B0_ksp_1;
            % initB0 = B0_ksp;
            initR2star = zeros(nSamples,nSamples,1);
        case 1
            % lastRAREdat = kSpace(:,:,:,:,:,1,i);
            % lastRAREtraj = traj(:,:,:,:,:,1);
            % bitmask2 = str2double(evalc("bart('bitmask 0 1')"));
            % initW = bart(sprintf('pics -S -R W:%i:0:0.0001 -i100 -e  -t',bitmask2),lastRAREtraj, lastRAREdat, sens);
            % initW = initW.*500*1000;
            initW = normalizeArray(lastSynthetic).*10;
            % initW = 2.*sqrt(initW.^2 +initW2.^2);
            initR2star = R2(:,:);
            initR2star = initR2star.*0.01.*2;
            initB0 = zeros(nSamples,nSamples,1);
        case 2
            initW = bart('extract 6 0 1',R_m0_1fe);
            % initW = initM0.*1;
%              initR2star = zeros(nSamples,nSamples,1);
            % initR2star = R2(:,:,i);
            initR2star = zeros(nSamples,nSamples,1);
            initB0 = initB0;
        case 3
            initW = bart('extract 6 0 1',R_m0_1fe);
            initR2star = zeros(nSamples,nSamples,1);
        case 4
            % initW = bart('extract 6 0 1',R_m0_1fe);
            % initW = initW.*1000/2;
            initW =ones(nSamples,nSamples).*10;
            initR2star = R2;
            initR2star = initR2star.*0.01.*1.3;
            initB0 = zeros(nSamples,nSamples,1);
        case 5
            initW = ones(nSamples,nSamples,1);
            initB0 = zeros(nSamples,nSamples,1);
            
            initR2star = R2(:,:,slice);
            initR2star = initR2star.*0.01;
    end
%     initR2star = R2(:,:,i);
    % scale by 0.01 so the inverse will be miliseconds, scale by 1.6 for
    % faster decay
%     initR2star = initR2star.*0.01.*2;
    init = bart('join 6',initW, initR2star, initB0);
    initMaps(:,:,:,:,:,:,:,slice) = init;
end
% come back to bart v09
functionDir = fileparts(mfilename('fullpath'));
bartPath = fullfile(functionDir,'../../../../tools/bart_v09/bart');
run(fullfile(bartPath, 'startup.m'));
end