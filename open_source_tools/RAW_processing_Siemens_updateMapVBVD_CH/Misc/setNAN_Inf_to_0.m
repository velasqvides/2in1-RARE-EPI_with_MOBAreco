function outputArray = setNAN_Inf_to_0( inputArray )

outputArray=inputArray;
outputArray(isnan(inputArray))=0;
outputArray(isinf(inputArray))=0; 
end
