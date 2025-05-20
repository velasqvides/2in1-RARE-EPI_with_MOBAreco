% function to check if input is string and converting it from hexadecimal
% to decimal
% this is needed, because reading the header of VB17 data with new mapVBVD
% version results in some entries containg hexadecimals
function y=hexa2deci(x)
    if ischar(x)
        y = hex2dec(x(3:end));
    else
        y=x;
    end
end