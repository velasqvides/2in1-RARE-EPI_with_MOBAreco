T2_test = readcfl('T2');
T2star_test = readcfl('T2star');
R2_test = readcfl('R2');
R2star_test = readcfl('R2star');
M0_test = readcfl('M0');


xRange = 25:245; 
yRange = 46:214;
T2_test_ = T2_test(xRange, yRange, :);
delayTime = 0.7;
fileName = 'T2_test.gif';
createGifForRelaxationMaps(T2_test_,175,'T2',delayTime,fileName,'T_2 (ms)');


xRange = 25:245; 
yRange = 46:214;
T2star_test_ = T2star_test(xRange, yRange, :);
fileName = 'T2star_test.gif';
createGifForRelaxationMaps(T2star_test_,125,'T1',delayTime,fileName,'T_2* (ms)');


R2prime_test = R2star_test - R2_test;
xRange = 25:245; 
yRange = 46:214;
R2prime_test_ = R2prime_test(xRange, yRange, :);
fileName = 'R2prime_test.gif';
createGifForRelaxationMaps(R2prime_test_,20,'T1',delayTime,fileName,'R_2'' (1/s)');

xRange = 25:245; 
yRange = 46:214;
M0_test_ = M0_test(xRange, yRange, :);
delayTime = 0.7;
fileName = 'M0_test.gif';
createGifForGrayImages(abs(M0_test_),delayTime,fileName,'PD');

