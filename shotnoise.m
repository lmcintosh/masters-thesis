function [process] = shotnoise(dt,T,N,figures)

%using mean = 0; std = 1.Sejnowski's 3 ms for tau. use dt = 0.5. T = 1000
%ms, resolution N = 100, just has to be small compared to length(u);
%adjusted dt accordingly.

%% Create filter

x = linspace(0,0.03,N);
filtt = x.*exp(-x*1000/3);
k = find(filtt>0.00001);
filtr = filtt(k);
xp = x(k);



%% Convolution

%desired_length = T/dt;

m = T/dt + length(filtr) - 1;
u = randn(m,1);

process = conv(u,filtr,'valid');
b = std(process);
process = process/b;

%actual_length = length(process);

if figures 
    
figure;
plot(process);
%figure; 
%plot(filtr)

end
