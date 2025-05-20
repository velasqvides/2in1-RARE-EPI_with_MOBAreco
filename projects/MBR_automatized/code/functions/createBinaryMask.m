function mask = createBinaryMask(sens, M0, recoSIZE)
bitmask = str2double(evalc("bart('bitmask 3')"));
sens = bart(sprintf('rss %i',bitmask),sens);
M0_rss_abs = abs(sens .* M0);
circleMask = createCircleMask([recoSIZE, recoSIZE], [recoSIZE + 1, recoSIZE + 1]/2, recoSIZE/2 - 3);
M0_rss_abs = bart(sprintf('resize -c 0 %i 1 %i',recoSIZE,recoSIZE),M0_rss_abs);
M0_rss_abs = M0_rss_abs .* circleMask;
mask = bart('threshold -B 0.4',M0_rss_abs);
se = strel('disk',100);
mask = imclose(mask,se);
end