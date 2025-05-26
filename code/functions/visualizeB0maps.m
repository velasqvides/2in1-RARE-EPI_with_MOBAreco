function visualizeB0maps(B0maps,minValue,maxValue)
    % B0maps(~isfinite(B0maps)) = 2000;
    nSlices = size(B0maps,3);
for j =1:nSlices
figure, imagesc(real((B0maps(:,:,j))),[minValue maxValue]);axis image; axis off; colormap('jet');
% Set NaN values to be black
        cmap = colormap;            % Get current colormap
        cmap(1, :) = [0, 0, 0];     % Set the first color to black
        colormap(cmap);             % Apply the modified colormap

        % Make NaN values appear black
        set(gca, 'Color', [0 0 0]);
end
end
