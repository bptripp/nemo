function testInput()     
    input = FunctionInput(@(t) t);

    tauRef = .0025;
    tauRC = .02;
    n = 500;
    intercepts = -1+2*rand(1,n);
    maxRates = 100+200*rand(1,n);
    V0 = 0;
    sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    test = CosinePopulation([1], sg, 'test');
    
    tau = .01;
    testInput = test.addTermination('input', tau, 1);
    testX = test.addOrigin('x', @(x) x);
    plotDecodedOrigin(test, 'x')
    
    n = Network(.001);
    n.addNode(input);
    n.addNode(test);
    n.addConnection(input.origins{1}, testInput);
    outputProbe = n.addProbe(testX, 'output');
    spikeProbe = n.addProbe(sg, SpikeGenerator.STATE_SPIKES);

    tic
    n.run(0, 1);
    toc
    plotProbe(outputProbe, .01)
    plotProbe(spikeProbe)
end