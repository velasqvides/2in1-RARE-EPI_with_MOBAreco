function [AEstimate, M] = lookupTableEstimate(signal,sigma,lookupTable)

% All images are scalled to sigma=1 to be able to use the same lookup table

M = abs(signal)/sigma;

Apoints = lookupTable.Apoints;
Mpoints = lookupTable.Mpoints;

if M >= max(Mpoints)
    AEstimate = sigma*M;
elseif M<= min(Mpoints)
    AEstimate = 0;
else

    larger = Mpoints >= M;
    firstLarger = find(larger,1);

    spacing = Apoints(2)-Apoints(1);

    AEstimate = Apoints(firstLarger-1) + spacing*(M - Mpoints(firstLarger-1))/(Mpoints(firstLarger) - Mpoints(firstLarger-1));

    AEstimate = sigma*AEstimate;
end

M = sigma*M;

end
