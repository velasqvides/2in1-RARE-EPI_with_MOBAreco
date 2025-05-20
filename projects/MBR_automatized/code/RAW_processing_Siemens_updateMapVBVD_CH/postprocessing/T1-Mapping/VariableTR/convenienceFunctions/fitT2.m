% author:   Jan Hentschel
% date:     22.03.2011

function[x,ssq,cnt] = fitT2(time, signal)
    %FITT2 uses LMFnlsq, startValues estimated
    sigNormalized=mat2gray(signal)-0.70;
    sigNormalizedAbs=abs(sigNormalized);
    [minValue,minIndex] = min(sigNormalizedAbs);
    sigStartingpoint=time(minIndex);
    
    
    
    startValues = [min(signal),max(signal)-min(signal),sigStartingpoint]; 
    residual = @(x) x(1)+x(2)*exp(-(time/x(3))) - signal;
    [x,ssq,cnt] = LMFnlsq(residual,startValues);
end