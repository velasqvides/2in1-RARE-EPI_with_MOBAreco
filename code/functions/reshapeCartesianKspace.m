function [kSpaceTEs] = reshapeCartesianKspace(kSpace)
 kSpaceTEs = permute(kSpace,[1 2 6 3 7 5 4]);
end
