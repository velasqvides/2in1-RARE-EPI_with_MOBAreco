function TE_T2 = createArrayTEs_RARE(protPara)
%UNTITLED24 Summary of this function goes here
%   Detailed explanation goes here
ETL_RARE = protPara.ETL_RARE;
ESL_RARE = protPara.ESP_RARE/1000;

echoTimes = (ESL_RARE:ESL_RARE:(ETL_RARE)*ESL_RARE)';
TE_T2 = ones(1,1,1,1,1,ETL_RARE);

TE_T2(1,1,1,1,1,:) = echoTimes; % first echo is discarded

end % end function