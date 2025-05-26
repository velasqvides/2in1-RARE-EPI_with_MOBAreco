function [kSpaceCorr] = removeK0drifts(kSpace,compon)
%REMOVEK0DRIFTS removes drifts from the kSpace center
%inputs: kspace - contains the lines in the first non-singleton dimension
%   compon - components to be removed, small positive whole number skalar or vector
% usage: [kSpaceCorr] = removeK0drifts(kSpaceEPI_16c_RV,1:2);

kSzOri=size(kSpace);
kSpace=squeeze(kSpace);
kSzSq=size(kSpace);
[u,sv,vt]=svd(squeeze(kSpace(kSzSq(1)/2,:,:)),0);
u=u-ones(kSzSq(2),1)*mean(u);
drift=u(:,compon)*sv(compon,compon)*vt(:,compon)';
kSpaceCorr = reshape(kSpace(:,:) - kSpace(:,:)/drift(:).'*drift(:).',kSzOri);
end