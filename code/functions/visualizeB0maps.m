function visualizeB0maps(B0maps,minValue,maxValue)
nSlices = size(B0maps,3);
for j =1:nSlices
    figure, imagesc(real((B0maps(:,:,j))),[minValue maxValue]);axis image; axis off; colormap('jet');
    cmap = colormap;            
    cmap(1, :) = [0, 0, 0];     
    colormap(cmap);             
    set(gca, 'Color', [0 0 0]);
end
end
