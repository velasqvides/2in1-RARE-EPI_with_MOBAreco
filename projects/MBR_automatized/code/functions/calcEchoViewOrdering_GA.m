function [echoViewOrdering] = calcEchoViewOrdering_GA(protPara)
% sort echoes to echo time according to view oredring scheme
% bit-rev, Theilmann et al., 2004
% in format: [no. of RARE echoes in ET  no. of EPI echoes in ET]
% ET = echo train


nRARE = protPara.ETL_RARE;
nEPI = protPara.ETL_EPI;
ETL = protPara.ETLfull;


if (protPara.ViewOrdering==3||protPara.ViewOrdering==4)
    if (nRARE == 4)
        echoViewOrdering = [1 5 3 7];
    elseif (nRARE == 8)
        echoViewOrdering = [1 9 5 13 3 11 7 15];
    elseif (nRARE == 16)
        echoViewOrdering = [1 17 9 25 5 21 13 29 3 19 11 27 7 23 15 31];
    elseif (nRARE == 32)
        echoViewOrdering = 1+2*[0 16 8 24 2 18 10 26 4 20 12 28 6 22 14 30 1 17 9 25 3 19 11 27 5 21 13 29 7 23 15 31];
    end
else
    echoViewOrdering = 1:2:2*nRARE;
end

%if (wipPara.interleaved==1) % not interleaved
echoViewOrdering = echoViewOrdering - floor(echoViewOrdering/2);
%end

% add EPI echo positions
% wipPara.interleaved = 1/2/3  -> "none"/"aligend"/"reversed"
if (protPara.interleaved==2||protPara.interleaved==3) % not interleaved
    if (protPara.interleaved==3)  % interleaved reversed
        %echoViewOrdering = [echoViewOrdering 2*nRARE:-2:2];
        echoViewOrdering = [echoViewOrdering 1:1:nEPI];
    else % interleaved aligned
        echoViewOrdering = [echoViewOrdering 2*(nRARE-1):2:4*(nRARE-1)];
    end
else % not interleaved
    echoViewOrdering = [echoViewOrdering (nRARE+1):ETL];
end
if (protPara.bSameSpokeInET==2)
    echoViewOrdering = 1:ETL;
end
end