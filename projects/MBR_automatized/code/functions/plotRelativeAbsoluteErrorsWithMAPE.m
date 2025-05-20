function mape= plotRelativeAbsoluteErrorsWithMAPE(reference,newMethod,binaryMask,maxError)
ref = reference.*binaryMask;
new = newMethod.*binaryMask;
errorMap = 100 .* abs( (ref - new) ./ ref );
errorMap(isnan(errorMap)) = 0;
errorMap(isinf(errorMap)) = 0; %20.37/13.46/8.81/28.22
figure,imagesc(errorMap,[0 maxError]); axis equal; axis off;
% theSum = sum(sum((errorMap(binaryMask>0))));
% N = nnz(theSum);
% mape=theSum/N;
mape = mean(errorMap(binaryMask>0));
end