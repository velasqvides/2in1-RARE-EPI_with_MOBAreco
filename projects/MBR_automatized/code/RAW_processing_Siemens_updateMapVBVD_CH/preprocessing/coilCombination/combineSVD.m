function [sensitivity,image] = combineSVD(data)

% [SENSITIVITY, IMAGE] = combineSVD(DATA)
%
% DATA is a Nx by Ny by Nz by Nreceivers by Nimages array of complex MRI data
%  All images and receivers should have the same noise, i.e. pre-whiten with noise covariance matrix
% SENSITIVITY is a Nx by Ny by Nz by Nreceivers complex sensitivity map array
% IMAGE is a Nx by Ny by Nz by Nimages complex image array

siz = size(data);
data = permute(data, [4 5 1 2 3]);
data = data(:,:,:);

sz = size(data);
image = zeros(sz(2),sz(3),'single');
sensitivity = zeros(sz(1),sz(3),'single');
for p = 1:sz(3)
   [U,S,V] = svd(data(:,:,p),'econ');

   % Image is first row of V**H, times largest eigenvalue
   %  Note that MATLAB svd returns V not V**H
   % Receiver sensitivity is first column of U
   % Force first image to be positive real, ie zero phase
   w = V(1)./abs(V(1));
   
   image(:,p) = conj(V(:,1)) .* S(1) .* w;
   sensitivity(:,p) = U(:,1) .* w;
end

image = permute(image, [2 1]);
image = reshape(image, siz([1 2 3 5]));

sensitivity = permute(sensitivity, [2 1]);
sensitivity = reshape(sensitivity, siz([1 2 3 4]));

