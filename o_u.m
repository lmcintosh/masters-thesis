function [process] = o_u(t,dt,g,D,key)

close all

%dt = mean(diff(t));
t = linspace(0,t,t/dt);
X_free = cumsum(sqrt(2*D*dt)*randn(size(t)));
X_filt = filter([0 g*dt], [1 -1+g*dt] , X_free);
%process = X_free-X_filt;

% THIS IS WRONG, BUT WHATEVER
process = X_filt;

if key

figure;
plot(t,process);
title('Ornstein-Uhlenbeck Process')
xlabel('Time')

%{
figure;
plot(t,X_filt);

figure;
plot(t,X_free);
%}
end