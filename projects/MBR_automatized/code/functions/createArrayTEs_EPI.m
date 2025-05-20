function TE_T2star = createArrayTEs_EPI(protPara)
%UNTITLED25 Summary of this function goes here
%   Detailed explanation goes here
ESL_EPI = protPara.ESP_EPI;
ETL_EPI = protPara.ETL_EPI;
firstEPI = protPara.ESP_EPIfirst;
rest_of_echoTimes = (ESL_EPI:ESL_EPI:(ETL_EPI)*ESL_EPI)';
rest_of_echoTimes = rest_of_echoTimes + firstEPI;
TE_T2star = zeros(1,1,1,1,1,ETL_EPI + 1);
TE_T2star(1,1,1,1,1,1) = 0;
TE_T2star(1,1,1,1,1,2) = firstEPI;
TE_T2star(1,1,1,1,1,3:end) = rest_of_echoTimes(1:end-1,:);
end