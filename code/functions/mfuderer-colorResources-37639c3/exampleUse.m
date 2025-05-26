% Load necessary matlab modules
addpath('/my/git/location');  


% Load the matlab file containing an example T1 map 
fn = 'sampleT1map.mat';
data = load(fn);
im = data.sampleT1;

loLev = 0.0;
upLev = 150.0;

% Call the relaxationColorMap function
[imClip, rgb_vec] = relaxationColorMap('T2', T2(:,:,1), loLev, upLev);

% Display the image using MATLAB
figure;
imshow(imClip, 'DisplayRange', [loLev, upLev], 'InitialMagnification', 'fit');
colormap(rgb_vec);
colorbar;