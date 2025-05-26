a=readcfl('T2star');
b=readcfl('T2star');

aa=a(:,:,[8 1 9 2 10 3 11 4 12 5 13 6 14 7]);
bb=b(:,:,[8 1 9 2 10 3 11 4 12 5 13 6 14 7]);

cc = zeros(256, 256, 28); % Preallocate the output array

% Interleave slices
for i = 1:14
    cc(:, :, 2*i-1) = aa(:, :, i); % Take the i-th slice from bb
    cc(:, :, 2*i) = bb(:, :, i);   % Take the i-th slice from aa
end

T2first=readcfl('T2');
T2starfirst=readcfl('T2star');
T2second=readcfl('T2');
T2starsecond=readcfl('T2star');

T2first_=T2first(:,:,[8 1 9 2 10 3 11 4 12 5 13 6 14 7]);
T2starfirst_=T2starfirst(:,:,[8 1 9 2 10 3 11 4 12 5 13 6 14 7]);
T2second_=T2second(:,:,[8 1 9 2 10 3 11 4 12 5 13 6 14 7]);
T2starsecond_=T2starsecond(:,:,[8 1 9 2 10 3 11 4 12 5 13 6 14 7]);

for i = 1:14
    T2final(:, :, 2*i-1) = T2first_(:, :, i); % Take the i-th slice from bb
    T2final(:, :, 2*i) = T2second_(:, :, i);   % Take the i-th slice from aa
end

for i = 1:14
    T2starfinal(:, :, 2*i-1) = T2starfirst_(:, :, i); % Take the i-th slice from bb
    T2starfinal(:, :, 2*i) = T2starsecond_(:, :, i);   % Take the i-th slice from aa
end

T2 = T2final(12:234,42:219,:);
T2star = T2starfinal(12:234,42:219,:);

config.dirToSave='/home/cstuser/Users/Velasquez/PhD_project/projects/MBR_automatized/reconstructions/volunteer_gaps3mm_noRRIshiming_301224/Merged_T2_T2star/T2';
figName  = 'T2';
finalDir = 'T2';
maptype  = 'T2'; 
maxValue = 175;
saveRelaxationMaps_newColorMaps(T2final, maxValue, figName, maptype, finalDir, config);

figName  = 'T2star';
finalDir = 'T2star';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps_newColorMaps(T2starfinal, maxValue, figName, maptype, finalDir, config);
