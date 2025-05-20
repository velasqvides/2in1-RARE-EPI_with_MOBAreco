function B1map=Interpolation_B1B0(B0map,B1dif,Curve,B0vec,B1vec)
% function B1map=Interpolation_B1B0(B0map,B1dif,Curve,B0vec,B1vec)
%
% This function takes B0 values from the B0 map and phase differences from
% the two B1+ sensitive acquisitions and looks up the interpolated value for
% relative B1+ (B1vec) based on the simulated B1+ sensitive phase difference.
%
% INPUT: B0Map = B0-Map
%        B1dif = B1 sensitive phase difference map
%        Curve = Interpolation curve based on bloch simulations
%        B0vec = range of B0 values simulated (vector)
%        B1vec = range of relative B1+ values simulated (vector)
%           
% MD interpolation simplified and about 3.5x faster 2013-06-06
% MD this function is still slow...need to vectorize it somehow
% 2013 Matthias Dieringer

B0vec=B0vec';
B1vec=B1vec';
B1map=zeros(size(B1dif));

for g=1:size(B1dif,1)
    for h=1:size(B1dif,2)
        if (~isnan(B0map(g,h)) && B1dif(g,h)~=0)
            %first: get the right row (B0)
            Atmp=interp2(B0vec,B1vec,Curve',B0map(g,h),B1vec, '*linear');
            %second: get right column (B1+) based on phase difference (Atmp)
            B1map(g,h)=interp1(Atmp,B1vec,B1dif(g,h),'linear');
        end
    end
end
end
