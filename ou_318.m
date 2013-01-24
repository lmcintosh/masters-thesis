function [process] = ou_318(naught,mu,tau,D,T,dt,noisetype,figures)

N = ceil(T/dt);

%% make stochF
stochF = zeros(N,1);

%% white noise
if (noisetype == 0)
 for n = 2:N;
  stochF(n) = randn;
 end
%% wiener process
elseif (noisetype == 1)
 for n = 2:N;
  stochF(n) = stochF(n-1) + randn;
 end 
end

%{
%% Simulate diffusion approximation
tseries = linspace(naught,T,N);
process = zeros(N,1);
noise1 = zeros(1,N-1);

stochF(1) = randn;
noise1(1) = (D*tau/2)*(1 - exp((-2*tseries(1))/tau))*stochF(1);

for n = 2:N;
    process(n-1) = mu + (exp(-tseries(n-1)/tau))*(naught - mu) + noise1(n-1);
    noise1(n) = (D*tau/2)*(1 - exp((-2*tseries(n))/tau))*stochF(n) + noise1(n-1);
end
process(N) = mu + (exp(-tseries(n-1)/tau))*(naught - mu) + noise1(N);
end
%}

%% Simulate stochastic differential equation
truth = zeros(N,1);
truth(1) = naught;
%mu = 0;

for n = 1:N-1;
    truth(n+1) = truth(n) - ((truth(n) - mu)/tau)*dt + (sqrt(D))*stochF(n);
end

process = truth;

if figures

figure;
plot(process)
title('Langevin')
xlabel('Time')

end