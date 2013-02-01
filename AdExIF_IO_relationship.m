% little script to match input/output relationship of AdExIF

time = 1000; % ms
delta = 0.5;
constants = linspace(0,1000);
currents = ones(time/delta,1)*constants; % we're running the simulation for 1000ms with delta = 0.5
adaptation = linspace(0,40);

rates = zeros(length(constants),length(adaptation));

for i = 1:length(constants)
    for j = 1:length(adaptation)
        [V w] = AdExIF(currents(:,i),adaptation(j));
        spiketrain = findSpikes(V,18);
        rates(i,j) = sum(spiketrain); % want firing rate as spikes/sec, and running simul. for 1 sec.
    end
end

figure; surf(adaptation,constants,rates) % counterintuitively, surf takes triplets (x(j), y(i), z(i,j))
ylabel('Current step magnitude')
xlabel('Adaptation parameter')
zlabel('Firing rate (Hz)')
title('Response of an Adaptive, Exponential Integrate-and-Fire Neuron to Current Steps')
