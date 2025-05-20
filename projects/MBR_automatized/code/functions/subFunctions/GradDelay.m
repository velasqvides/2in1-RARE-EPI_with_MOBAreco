 %it work perfectly by using the same spoke and then adding a delay to it!

close all
spokeA = gdelcalib.gdelcalib_0(:,1,1,1,1);
spokeB = gdelcalib.gdelcalib_180(:,1,1,1,1);
fft_spokeA = fftshift(fft((spokeA)));
fft_spokeB = fftshift(fft(flipud((spokeB))));
g = fft_spokeA.*conj(fft_spokeB);
slope = angle(g);
figure, plot(slope)
cropedSig = slope(189:325);
p = polyfit(1:length(cropedSig),cropedSig,1);
GradDelayCal = -p(1)*length(spokeA)/(4*pi)


g1 = circshift(g,-1);
scenario2 = g(189:325).*g1(189:325);
slope1 = angle(sum(scenario2));
GradDelayCal = -slope1*length(spokeA)/(4*pi)
figure, plot(slope1)



%% another try. It worked using same spoke
spokeA = gdelcalib.gdelcalib_0(:,1,1,1,1);
delay = 0.3;
N = length(spokeA);
k = fftshift([0:floor(N/2)-1 floor(-N/2):-1])'; 
phasefactor = exp(-1i * 2 * pi * delay *  k / N);

fft_spokeA = fftshift(ifft(ifftshift(spokeA)));
spokeB = spokeA;
fft_spokeB = fftshift(ifft(ifftshift(spokeB)));
fft_spokeB_linearPhase = fft_spokeB.*phasefactor;
spokeB_moved = fftshift(fft(ifftshift(fft_spokeB_linearPhase)));
fft_spokeB_moved = fftshift(ifft(ifftshift(spokeB_moved)));
g = fft_spokeA.*conj(fft_spokeB_moved);
slope = angle(g);
cropedSig = slope(189:325);
p = polyfit(1:length(cropedSig),cropedSig,1);
GradDelayCal = -p(1)*length(spokeA)/(2*pi)






