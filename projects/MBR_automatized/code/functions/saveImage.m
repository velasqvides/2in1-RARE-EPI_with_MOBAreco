function saveImage(T2maps,maxValue,figName,config)
% at the moment T2 MOBA works for bartv07 only
dirToSave = config.dirToSave;
finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave,'images');
nSlices = size(T2maps,3);
for j =1:nSlices
figure, imagesc((T2maps(:,:,j)),[0 maxValue]);axis equal; axis off; 
set(gcf, 'Units', 'pixels', 'Position', [0, 0, size(T2maps, 2), size(T2maps, 1)]);
% Set the axes position to fill the entire figure
set(gca, 'Position', [0, 0, 1, 1]);
sliceN = strcat(figName, '_', num2str(j));
print(fullfile(filePath,sliceN),'-dpng','-r0');
end
end

