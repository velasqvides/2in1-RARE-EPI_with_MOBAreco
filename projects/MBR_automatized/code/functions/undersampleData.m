function [kSpace_6, traj_6] = undersampleData(kSpace, traj, protPara, ETL)

%% 1. Give the acceleration factor retrieve some reconstruction parameters
accFactor = 6;
nSpokes = protPara.nSpokes;
readoutSamples = protPara.baseRes*protPara.oversamplingFactor;
spokesPerKspace = ceil(nSpokes/accFactor);
nSpokesReco = spokesPerKspace*ETL;
initialIndex = 1:accFactor:nSpokes;
nChannels = protPara.nChannels;
nSlices = protPara.nSlices;

%% 3. Name the undersmapled data

traj_6 = zeros(3,readoutSamples,spokesPerKspace,1,1,ETL);
kSpace_6 = zeros(1,readoutSamples,spokesPerKspace,nChannels,1,ETL,nSlices);

%% 4. retrospectively undersample the data
toSum = [0 3 1 5 2 4 0 3 1 5 2 4 0 3 1 5 2 4 0 3 1 5 2 4]; % for acc=6;
indexOrder = [];
for i=1:length(toSum)
    indexOrder = [indexOrder ( initialIndex + toSum(i) )];
    indexOrder(indexOrder>nSpokes) = [];
end
finalIndexOrder = indexOrder(1:nSpokesReco);
size(finalIndexOrder);

sindx = 1; eindx = spokesPerKspace;

for i=1:ETL
    traj_6(:,:,:,:,:,i) = traj(:,:,finalIndexOrder(sindx:eindx),:,:,i);
    sindx = sindx + spokesPerKspace; eindx = eindx + spokesPerKspace;
end


for i=1:nSlices
    sindx = 1; eindx = spokesPerKspace;
for j=1:ETL
    kSpace_6(:,:,:,:,:,j,i) = kSpace(:,:,finalIndexOrder(sindx:eindx),:,:,j,i);
    sindx = sindx + spokesPerKspace; eindx = eindx + spokesPerKspace;
end
end