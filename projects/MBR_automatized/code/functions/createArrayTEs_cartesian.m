function TE_T2 = createArrayTEs_cartesian(protPara)
%UNTITLED24 Summary of this function goes here
%   Detailed explanation goes here
ETL = protPara.ETL;
ESL = protPara.ESP/1000;

echoTimes = (ESL:ESL:(ETL)*ESL)';
TE_T2 = ones(1,1,1,1,1,ETL);

TE_T2(1,1,1,1,1,:) = echoTimes; % first echo is discarded

end % end function