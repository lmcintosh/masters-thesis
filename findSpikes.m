function [spikeTrain] = findSpikes(V,threshold)
% INPUTS: V is raw voltage trace, and a spike is called when V >= threshold.
% OUTPUS: 0 or 1 spiketrain in the same dimensions as V.

spikeTrain = (V >= threshold);
