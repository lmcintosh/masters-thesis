function [process, trueVariance] = wiener(startingvalue, drift, variance, N, fig)
% INPUTS: startingValue, drift, variance, N, fig
% OUTPUT: process, trueVariance
% The Wiener process X_t = drift*t + sqrt(variance)*W_t has variance drift^2*t^2 - drift*t + variance^2*t

truth = zeros(N,1);
truth(1) = startingvalue;

for n = 1:N-1;
    truth(n+1) = truth(n) + drift + sqrt(variance) * randn;
end
process = truth;

trueVariance = (drift^2)*(N^2) - (drift*N) + variance*N;
% Note: inaccurate for small N

if (fig)
    figure;plot(process)
end
