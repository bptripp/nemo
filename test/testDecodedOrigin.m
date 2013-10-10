% Tests for DecodedOrigin

% merge ... 
tauRef = .0025;
tauRC = .02;
n = 10;
intercepts = -1+2*rand(1,n);
maxRates = 100+200*rand(1,n);
V0 = 0;
sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
cp = CosinePopulation(1, sg, 'test');
do = addOrigin(cp, 'output', @(x) [x; x.^2]);

addMerge(do, [1 2 3]);
assert(sum(do.decoders(:,11) - sum(do.decoders(:,1:3), 2)) == 0)
removeMerge(do, 11);
assert(size(do.decoders, 2) == 10)

disp('Test OK')
