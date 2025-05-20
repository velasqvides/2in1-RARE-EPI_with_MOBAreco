function diff=phaseDiff(z1,z2)
%% Olli Kraus (19.6.2013)
% function uses atan2
%
% input: z1, z2 must be complex arrays of same size (can be 2D or more!)
% output: diff is a real array of same size as z1 or z2 with values [-pi +pi]
%
if isreal(z1) && isreal(z2)
    disp('input must be complex!')
    diff=[];
    return
else
diff=atan2(imag(z1 .* conj(z2)), real(z1 .* conj(z2)));
end