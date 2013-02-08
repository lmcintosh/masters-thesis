function [process] = wiener(startingvalue, sigma, N, fig)
% INPUTS: startingValue, sigma, N, fig
% OUTPUT: process

truth = zeros(N,1);
truth(1) = startingvalue;

for n = 1:N-1;
    truth(n+1) = truth(n) + sigma * randn;
end
process = truth;

if (fig)
    figure;plot(process)
end
