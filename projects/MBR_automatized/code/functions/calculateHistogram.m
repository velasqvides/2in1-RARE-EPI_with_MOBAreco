map=T2_carl;
bin_edges = 40:2:200;
histogram_values = histcounts(map, bin_edges);
figure;
bar(bin_edges(1:end-1), histogram_values);
title('Histogram of Image T2_carl');
xlabel('Pixel Value');
ylabel('Frequency');