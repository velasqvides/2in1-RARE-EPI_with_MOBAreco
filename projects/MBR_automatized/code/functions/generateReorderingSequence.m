function order = generateReorderingSequence(n)
order = zeros(1, n);
if mod(n, 2) == 0
    order(1:2:end) = 1:(n/2);
    order(2:2:end) = (n/2 + 1):n;
else
    order(1:2:end) = 1:ceil(n/2);
    order(2:2:end) = (ceil(n/2) + 1):n;
end
end