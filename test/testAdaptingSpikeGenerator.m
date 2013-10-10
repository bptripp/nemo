% TODO: use unit testing framework
function testAdaptingSpikeGenerator() 
    input = FunctionInput(@(t) 1);

    tauRef = .0025;
    tauRC = .02;
    n = 50;
    intercepts = -1+2*rand(1,n);
    maxRates = 100+200*rand(1,n);
    V0 = 0;
    tauAdapt = .1; % uniform
    Ginc = 1 + rand(1,n);
    Gad0 = 0;
    sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0, tauAdapt, Ginc, Gad0);
    test = CosinePopulation([1], sg, 'test');
    
    tau = .01;
    testInput = test.addTermination('input', tau, 1);
    
%     dt = .001;
%     T = 1;
%     time = 0:dt:T;
% 
%     spikes = zeros(n, length(time)-1);
%     for j = 1:(length(time)-1)
%         spikes(:,j) = sg.run(1, time(j), time(j+1), 1);
%     end        

    n = Network(.001);
    n.addNode(input);
    n.addNode(test);
    n.addProjection(input.origins{1}, testInput);
    spikeProbe = n.addProbe(sg, SpikeGenerator.STATE_SPIKES);

    tic
    n.run(0, 1);
    toc
    plotProbe(spikeProbe)
end
