function [ w ] = tukey( width, alpha )
%TUKEY creates a Tukey window of a given width
%   
%   USAGE:
%       [ w ] = tukey( width, alpha )
%
%       width = window width
%       alpha = parameter that defines part of the window center that is 1
%               for alpha = 1, window becomes a Hann window
%               for alpha = 0, window becomes rectangular
%               default = 1/3
%           
% for details see https://en.wikipedia.org/wiki/Window_function#Tukey_window

% 2015-11-12, Till Huelnhagen, created

    if nargin<2
        alpha = 1/3;
    end
    
    w = zeros(width,1);

    for n = 0:width-1
        if (n >= 0 && n <= (alpha*(width-1) / 2))
            w(n+1) = 0.5 * (1 + cos(pi*(2*n / (alpha*(width-1)) - 1)));
        elseif (n >= (alpha*(width-1) / 2) && n <= ((width-1)*(1-alpha/2)))
            w(n+1) = 1;
        else
            w(n+1) = 0.5 * (1 + cos(pi*(2*n / (alpha*(width-1)) - 2/alpha + 1)));
        end
    end
    w(isnan(w))=1;
end
