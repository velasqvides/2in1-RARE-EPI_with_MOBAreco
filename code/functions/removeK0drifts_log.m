function [kSpaceCorr] = removeK0drifts_log(kSpace,compon)
%REMOVEK0DRIFTS removes drifts from the kSpace center
%inputs: kspace - contains the lines in the first non-singleton dimension
%   compon - components to be removed, small positive whole number skalar or vector
% usage: [kSpaceCorr] = removeK0drifts_log(kSpaceEPI_16c_RV,2);

kSzOri=size(kSpace);
kSpace=squeeze(kSpace);
kSzSq=size(kSpace);
[u,sv,vt]=svd([log(abs(squeeze(kSpace(kSzSq(1)/2,:,:)))) angle(squeeze(kSpace(kSzSq(1)/2,:,:)))],0);
%u(:,1)=u(:,1)-mean(real(u(:,1))); %ones(kSzSq(2),1)*
sv=diag(sv);
vt = vt(1:end/2,1:length(sv)) + 1i*vt((end/2+1):end,1:length(sv));
drift=double(u(:,compon)*diag(sv(compon))*vt(:,compon).');

%driftFitR=real(log(double(kSpace(:,:))))/real(drift(:)');
%driftFitI=imag(log(double(kSpace(:,:))))/imag(drift(:).');
driftFitR=real(log(double(kSpace(:,:)))/drift(:).');
%driftFitR=real(driftFit);
for ik=reshape([(kSzSq(1)/2+1):(kSzSq(1)-1); (kSzSq(1)/2-1):-1:1],1,[]) % starting in the center
    driftFitR(ik)= fminbnd(@(x) fitPhase(x,log(double(kSpace(ik,:))),drift(:).'),-1,1);
%     disp([ik/kSzSq(1) real(driftFit(ik)) driftFitR(ik)])
end
kSpaceCorr = single(exp(reshape(log(double(kSpace(:,:))) - driftFitR*drift(:).',kSzOri)));

function [devia] = fitPhase(driftFit,logKSpace,drift)
	devia=logKSpace-driftFit*drift;
	devia=real(devia)+1i*(rem(imag(devia)+pi,2*pi)-pi);
	devia=sum(sum(abs(devia).^2));
end
end