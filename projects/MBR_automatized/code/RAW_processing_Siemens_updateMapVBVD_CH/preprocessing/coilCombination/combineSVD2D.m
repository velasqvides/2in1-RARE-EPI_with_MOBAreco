function [sensitivity,image] = combineSVD2D(data)
% [SENSITIVITY, IMAGE] = combineSVD(DATA)
% Calculates coil sensitivities and combines receiver channels of MR images
%
% DATA is a Nx by Ny by Nreceivers by Nimages array of complex MRI data
% All images and receivers should have the same noise, i.e. pre-whiten with
% noise covariance matrix.
% SENSITIVITY is a Nx by Ny by Nreceivers complex sensitivity map array
% IMAGE is a Nx by Ny by Nimages complex image array
%
% Code by Martyn Klassen 2013 (Proc Intl. Soc. Magn. Reson. Imag., p.3739)
% adapted by Matthias Dieringer 2013 (2D Version)

siz = size(data);
data = permute(data, [3 4 1 2]);
data = data(:,:,:);

sz = size(data);
image = zeros(sz(2),sz(3));
sensitivity = zeros(sz(1),sz(3));
for p = 1:sz(3)
   [U,S,V] = svd(data(:,:,p),'econ');

   % Image is first row of V**H, times largest eigenvalue
   %  Note that MATLAB svd returns V not V**H
   % Receiver sensitivity is first column of U
   % Force first image to be positive real, ie zero phase
   w = V(1)./abs(V(1));
   
   image(:,p) = conj(V(:,1)) .* S(1) .* w;
   sensitivity(:,p) = U(:,1) .* w;
   
%       w = V(1)./abs(V(1));
%    
%    image(:,p) = conj(V(:,1)) .* S(1) .* V(1);
%    sensitivity(:,p) = U(:,1) .* V(1);
end

image = permute(image, [2 1]);
image = reshape(image, siz([1 2 4]));

sensitivity = permute(sensitivity, [2 1]);
sensitivity = reshape(sensitivity, siz([1 2 3]));

