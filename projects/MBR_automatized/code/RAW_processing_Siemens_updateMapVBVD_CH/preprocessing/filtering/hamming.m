% HAMMING creates a hamming windowing function of length L
%
%   USAGE:
%       hamming(L)
%

function w = hamming(L)
    L=L+1;
    w=(0.54-0.46*cos(2*pi*(1:L)/L))';
    w=w(1:L-1);
end