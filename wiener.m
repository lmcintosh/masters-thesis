function [process] = wiener(startingValue, drift, trueVariance, N, fig)
% INPUTS: startingValue, drift, variance, N, fig
% OUTPUT: process, trueVariance
% The Wiener process X_t = drift*t + sqrt(variance)*W_t has variance drift^2*t^2 - drift*t + variance^2*t

%trueVariance = (drift^2)*(N^2) - (drift*N) + variance*N;
% Note: inaccurate for small N
sigma = sqrt((trueVariance - (drift^2)*(N^2) + (drift*N))/N);

displacements = drift + sigma*randn(N,1);
process = startingValue + cumsum(displacements);


if (fig)
    figure;plot(process)
end
