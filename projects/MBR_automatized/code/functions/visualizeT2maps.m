function visualizeT2maps(T2maps,maxValue)
nSlices = size(T2maps,3);
for j =1:nSlices
figure, imagesc((T2maps(:,:,j)),[0 maxValue]);axis image; axis off; 
end
end