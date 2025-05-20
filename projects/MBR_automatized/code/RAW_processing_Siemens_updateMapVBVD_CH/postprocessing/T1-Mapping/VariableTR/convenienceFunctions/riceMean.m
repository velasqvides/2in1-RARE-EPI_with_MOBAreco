function mean = riceMean(A,sigma)

    x = -A.^2/(2*sigma^2);

    mean = sigma*sqrt(pi/2)*exp(x/2).*((1-x).*besseli(0,-x/2) - x.*besseli(1,-x/2));

end