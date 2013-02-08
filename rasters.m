function [] = rasters(spikes,T)
% INPUTS: spikes, T
% OUTPUTS: none

numbers = size(spikes,1);
ts = 0:numbers-1;
M = size(spikes,2);

if nargin~=2
    T = numbers;
end

k = 0;
figure;
axis([0,T,0,M])
%colors = varycolor(M);
for j = 1:M;
    spiketrain = ts(find(spikes(:,j)));
    if isempty(spiketrain)
        k = k+1;
    else
        for i = 1:length(spiketrain)
            line([spiketrain(i),spiketrain(i)],[k,k+1],'Color','k');
            hold on
        end
        k = k+1;
    end
    % title(['Neuron ',num2str(j)]);
end
ylabel('Neuron #')
title('Raster Plot')