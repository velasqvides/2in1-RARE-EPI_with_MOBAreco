function [initMaps] = createInitMpasForT2star_m3_v09(kSpace, traj, TEs, sensRARE, R2, synthesizedRAREimages, binaryMaskRARE, whichInit)
% B0_ksp = bart(sprintf('fft -u %i', bitmask2), initB0_fessler);
% only needt o get an idea of current site for moba for R2*mappig
TE_1 = TEs(:,:,:,:,:,1:3);
traj_1 = traj(:,:,:,:,:,1:3,1);
kdat_1 = kSpace(:,:,:,:,:,1:3,1);
dummy = bart(sprintf('moba -g -G -m0 -i1 -o1  -t'),traj_1, kdat_1, TE_1);
dummy = bart('extract 6 2 3', dummy);
%
nSamplesRows = size(dummy, 1);
nSamplesCols = size(dummy, 2);
nSlices = size(kSpace,7);
initMaps = zeros(nSamplesRows,nSamplesCols,1,1,1,1,3,nSlices);
% at the moment T2 MOBA works for bartv07 only
changeBartVersion(9)
for slice = 1:nSlices
    switch whichInit
        case 0
            initW = ones(nSamplesRows,nSamplesCols,1);    
            initR2star = zeros(nSamplesRows,nSamplesCols,1);
            initB0 = zeros(nSamplesRows,nSamplesCols,1);
        case 1
            toPadRows = (nSamplesRows - size(R2, 1))/2;
            toPadCols = (nSamplesCols - size(R2, 2))/2;
            lastSynthetic = synthesizedRAREimages(:,:,:,:,:,end,slice).*binaryMaskRARE(:,:,slice);
            initW = normalizeArray(lastSynthetic).*10;
            initW = padarray(initW,[toPadRows toPadCols]);
            initR2star = R2(:,:,slice).*binaryMaskRARE(:,:,slice);
            initR2star = initR2star.*0.01.*2;
            % initW = padarray(initR2star,[toPadRows toPadCols]);
            % initR2star = zeros(nSamples,nSamples,1);
            initR2star = padarray(initR2star,[toPadRows toPadCols]);
            initB0 = zeros(nSamplesRows,nSamplesCols,1);
        case 2
            TE_1 = TEs(:,:,:,:,:,1:3);
            traj_1 = traj(:,:,:,:,:,1:3,slice);
            kdat_1 = kSpace(:,:,:,:,:,1:3,slice);
            init_m0 = bart(sprintf('moba -g -O -G -m0 -i6 -o%f --fat_spec_0 -R3 -t',1),traj_1, kdat_1, TE_1);
            initB0 = bart('extract 6 2 3', init_m0);
            initW = bart('extract 6 0 1',init_m0);
            initR2star = zeros(size(initW,1),size(initW,2),1);
        case 3
            toPadRows = (nSamplesRows - size(R2, 1))/2;
            toPadCols = (nSamplesCols - size(R2, 2))/2;
            lastSynthetic = synthesizedRAREimages(:,:,:,:,:,end,slice).*binaryMaskRARE(:,:,slice);
            initW = normalizeArray(lastSynthetic).*10;
            initW = padarray(initW,[toPadRows toPadCols]);
            initR2star = R2(:,:,slice).*binaryMaskRARE(:,:,slice);
            initR2star = initR2star.*0.01.*2;
            initR2star = padarray(initR2star,[toPadRows toPadCols]);
            TE_1 = TEs(:,:,:,:,:,1:3);
            traj_1 = traj(:,:,:,:,:,1:3,slice);
            kdat_1 = kSpace(:,:,:,:,:,1:3,slice);
            init_m0 = bart(sprintf('moba -g -O -G -m0 -i6 -o%f --fat_spec_0 -R3 -t',1),traj_1, kdat_1, TE_1);
            initB0 = bart('extract 6 2 3', init_m0);
    end
    init = bart('join 6',initW, initR2star, initB0);
    initMaps(:,:,:,:,:,:,:,slice) = init;
end
% come back to bart v09
changeBartVersion(9)
end