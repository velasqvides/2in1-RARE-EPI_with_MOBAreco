function full_kspace = Grappa2D(RAW,ref,R)
% function full_kspace = Grappa2D(RAW,ref,R,bs)
% This is 2D GRAPPA reconstruction
%
% INPUT:
% RAW       : subsampled RAW data (Array of [Freq x Phase x Coils]) with
%             zero-lines(!)
% ref       : fully sampled reference scan
% R         : acceleration (reduction) factor
% bs        : block size for calibration, default is 3 (must be odd)
%
% OUTPUT    : Fully reconstructed k-space
%
% Todo:  -There is still a bug with block sizes other than 3. You could
%         use 5 blocks e.g., but the number of lines in the reconstructed
%         k-space are then not complete
%        -this code is already quite fast, but loops can be further
%         vectorized to speed up
%
% 06-2013 Matthias Dieringer

% Perform some checks
% default bs is 3
% if nargin<4
    bs=3;
% end

% if no acceleration factor is given, try to find out
if nargin<3
    R=diff(find(RAW(1,:,1),2));
end

% reference data size
[NXref,NYref,nc]=size(ref);
% RAW data size
[NXraw,NYraw,~]=size(RAW);
% downwards rounded block size / 2
bs2=fix(bs/2);

% 1 First calculate GRAPPA weights
src=zeros((NYref-(bs-1)*R)*(NXref-(bs-1)),nc*bs*bs);
targ=zeros((NYref-(bs-1)*R)*(NXref-(bs-1)),nc*(R-1));

cnt=0;                          
for yind=1:NYref-(bs-1)*R
    for xind=bs2+1:NXref-bs2
        cnt=cnt+1;
        src(cnt,:)=reshape(ref(xind-bs2:xind+bs2, yind:R:yind+(bs-1)*R, :), 1, nc*bs*bs);
        targ(cnt,:)=reshape(ref(xind, yind+(bs2-1)*R+1 : yind+(bs2)*R-1, :), 1, nc*(R-1));
    end
end
ws=src\targ;

% 2 EXTEND MATRIX
% Extended matrix (readout direction) depending on kernel size
% this is needed to take care of the edges
nxExt = NXraw+2*bs2;

% Extended matrix (phase encoding direction) depending on kernel size
% this is needed to take care of the edges
% find first non-zero line to determine the number of lines to add
AddLinesFront = R-find(RAW(1,1:R,1),1)+1;
% calculate the number of lines to add to the back
AddLinesBack  = mod( AddLinesFront+NYraw, R) + (bs2-1)*R + 1;
% finally add those lines
nyExt = AddLinesFront+AddLinesBack+NYraw;

% 3 GRAPPA RECONSTRUCTION
% allocate memory
full_kspace=zeros(nxExt, nyExt, nc);
% copy raw data into enlarged array
full_kspace((bs2+1):(end-bs2), AddLinesFront+1:end-AddLinesBack, :) = RAW;
% perform fitting
for yind = 1:R:nyExt-(bs-1)*R
    for xind = bs2+1 : nxExt-bs2
        src=reshape(full_kspace(xind-bs2:xind+bs2, yind:R:yind+(bs-1)*R, :), 1, nc*bs*bs);
        full_kspace(xind,yind+(bs2-1)*R+1:yind+bs2*R-1,:) = reshape(src*ws,[R-1 nc]);
    end
end

% finally extract fully reconstructed k-space without leading and trailing
% zeros
full_kspace(:,~any(full_kspace,1))=[];
full_kspace(~any(full_kspace,2),:)=[];
full_kspace=reshape(full_kspace,[NXraw,NYraw,nc]);