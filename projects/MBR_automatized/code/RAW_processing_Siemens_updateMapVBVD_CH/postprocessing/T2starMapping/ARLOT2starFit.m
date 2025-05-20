% Fast T2* fitting approach based on ARLO
% (Pei et al., Proc. Intl. Soc. Mag. Reson. Med. 22 (2014), #0339)
%
% USAGE:
%
%   [t2starMap] = ARLOT2starFit(mag, deltaTE)
%
%   mag = array of magnitude images
%   deltaTE = inter echo spacing
%   ECHO DIMENSION MUST BE LAST IN ARRAY
%   deltaTE must be given in seconds

% 2014-05-29, Till Huelnhagen

function [t2starMap] = ARLOT2starFit(mag, deltaTE)

% reshape data to have echo position first
permuteVec=1:ndims(mag);
permuteVec=circshift(permuteVec,[0 1]);
mag=permute(mag,permuteVec);

nEcho = size(mag,1);

siz=size(mag);
sum1=zeros(siz(1,2:end));
sum2=sum1;
% calculate sums
for eco = 1:nEcho - 2
    sum1(:) = sum1(:) + ((deltaTE/3 * (mag(eco,:) + 4*mag(eco+1,:) + mag(eco+2,:))).^2)';
    sum2(:) = sum2(:) + ((deltaTE/3 * (mag(eco,:) + 4*mag(eco+1,:) + mag(eco+2,:))).*(mag(eco,:)-mag(eco+2,:)))';
end
% apply fitting function
t2starMap = (sum1 + deltaTE/3 * sum2) ./ (deltaTE/3 * sum1 + sum2);

end
