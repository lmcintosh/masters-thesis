function [V w] = AdExIF(current,T)
% current I is matrix with size T/delta x ensemble size and T is the stimulus duration in seconds
N = size(current,1);
M = size(current,2);

%% Neuron Parameters
C = 281; % capacitance in pF ... this is 281*10^(-12) F
g_L = 30; % leak conductance in nS
E_L = -70.6; % leak reversal potential in mV ... this is -0.0706 V
delta_T = 2; % slope factor in mV
V_T = -50.4; % spike threshold in mV
tau_w = 144; % adaptation time constant in ms
V_peak = 20; % when to call action potential in mV
b = 0.0805; % spike-triggered adaptation

loadPhysicalConstants;

%% Init variables
V = zeros(N,M);
w = zeros(N,M);

%% Boltzmann distributed initial conditions
sigma_sq = 1/(beta*(C*10^(-12))); % make sure C is in F
% This variance is on the order of 1.5*10^(-11)
V(1,:) = E_L/1000 + sqrt(sigma_sq)*randn(1,M); % make sure E_L is in V
V(1,:) = V(1,:)*1000; % convert initial conditions into mV
%V(1,:) = -70 + (-40+70).*rand(1,M);
w(1,:) = 0; % initial condition is w = 0

%% Loop over the 1 <= n <= N time intervals
for n = 1:N-1
    V(n+1,:) = V(n,:) + (delta/C)*(spiking(V(n,:),g_L,E_L,delta_T,V_T)-w(n,:)+ I(n+1,:));
    w(n+1,:) = w(n,:) + (delta/tau_w)*(a*(V(n,:)-E_L) - w(n,:));
    
    % spiking mechanism
    alreadySpiked = (V(n,:) == V_peak);
    V(n+1,alreadySpiked) = E_L;
    w(n+1,alreadySpiked) = w(n,alreadySpiked) + b;

    justSpiked = (V(n+1,:) > V_peak);
    V(n+1,justSpiked) = V_peak;
end
