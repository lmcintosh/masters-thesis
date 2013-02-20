function [process] = wiener(startingValue, drift, sigma, N, T, fig)
% INPUTS: startingValue, drift, variance, N, T, fig
% OUTPUT: process, trueVariance
% The Wiener process X_t = drift*t + sqrt(variance)*W_t has variance drift^2*t^2 - drift*t + variance^2*t

%trueVariance = (drift^2)*(N^2) - (drift*N) + variance*N;
% Note: inaccurate for small N
dt = T/N;

displacements = drift + sigma*sqrt(dt)*randn(N,1);
process = startingValue + cumsum(displacements);


if (fig)
    figure;plot(process)
end
