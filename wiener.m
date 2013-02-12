function [process] = wiener(startingvalue, drift, trueVariance, N, fig)
% INPUTS: startingValue, drift, variance, N, fig
% OUTPUT: process, trueVariance
% The Wiener process X_t = drift*t + sqrt(variance)*W_t has variance drift^2*t^2 - drift*t + variance^2*t

truth = zeros(N,1);
truth(1) = startingvalue;

%trueVariance = (drift^2)*(N^2) - (drift*N) + variance*N;
% Note: inaccurate for small N
sigma = sqrt((trueVariance - (drift^2)*(N^2) + (drift*N))/N);

for n = 1:N-1;
    truth(n+1) = truth(n) + drift + sigma * randn;
end
process = truth;

if (fig)
    figure;plot(process)
end
