% Test for merging of LIF spike generators for reduced simulations.

tauRef = .0025;
tauRC = .02;
n = 10;
intercepts = -1+2*rand(1,n);
maxRates = 100+200*rand(1,n);
V0 = 0;
sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);

ind = [1 2];
addMerge(sg, ind);
assert(length(sg.scales) == n+1);
assert(sg.mergedCount(n+1) == 2)

drive = -1:.05:1;
ratesInd = getRates(sg, repmat(drive, length(ind), 1), 0, 0, ind);
ratesMerge = getRates(sg, drive, 0, 0, n+1);
% figure, hold on, plot(drive, ratesInd, 'k'), plot(drive, ratesMerge, 'r.')

dt = .001;
time = dt:dt:.5;

spikes = zeros(3, length(time));
for i = 1:length(time)
    spikes(:,i) = integrate(sg, repmat(drive(end), 3, 1), time(i)-dt, time(i), [1; 2; 11]);
end
assert( abs((spikes(1)+spikes(2))/2 - spikes(3)) <= 1, 'Merged neuron not spiking at max rate that is average of original neurons')
% sum(spikes, 2)
% plot(time, spikes + [0; 1; 2] * ones(size(time)))

removeMerge(sg, 11);
assert(length(sg.scales) == n);

disp('Test OK')
