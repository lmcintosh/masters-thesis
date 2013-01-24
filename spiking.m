function [ f_V ] = spiking(V,g_L,E_L,delta_T,V_T)
%SPIKING f(v), the spiking behavior of the neuron
%   nonlinear spiking properties as a function of voltage

M = length(V);
f_V = zeros(1,M);

for i = 1:M;
    f_V(i) = -g_L*(V(i) - E_L) + g_L*delta_T*exp((V(i) - V_T)/delta_T);
end
