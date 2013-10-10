function testIntegrator()     
    input = FunctionInput(@(t) 1);

    tauRef = .0025;
    tauRC = .02;
    n = 500;
    intercepts = -1+2*rand(1,n);
    maxRates = 100+200*rand(1,n);
    V0 = 0;
    sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    integrator = CosinePopulation([1], sg, 'integrator');
    
    tau = .1;
    integratorInput = integrator.addTermination('input', tau, tau);
    integratorFeedback = integrator.addTermination('feedback', tau, 1);
    integratorX = integrator.addOrigin('x', @(x) x);
    
    n = Network(.001);
    n.addNode(input);
    n.addNode(integrator);
    n.addConnection(input.getOrigin(), integratorInput);
    n.addConnection(integratorX, integratorFeedback);
    
    n.addProbe(input.origins{1}, 'output');
    n.addProbe(integratorX, 'output');

    n.runAndPlot(0, 1);
end