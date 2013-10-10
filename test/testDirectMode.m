% tests: 
% - subnetwork OK 
% - arbitrary node
% - multiple terminations OK
% - feedback OK
% - unused origin OK
% - functioninput OK
% - several interconnected populations OK
% - multidimensional populations OK
% - multidimensional terminations, origins OK
function testDirectMode()
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
    
    inputProbe = n.addProbe(input.origins{1}, 'output');
    outputProbe = n.addProbe(integratorX, 'output');
    spikeProbe = n.addProbe(sg, SpikeGenerator.STATE_SPIKES);
    
    % add another population that calculates a product and test again ... 
    sg2 = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    product = CosinePopulation([1; 1], sg2, 'product');
    productA = product.addTermination('input', tau, [1; 0]);
    productB = product.addTermination('integrated', tau, [0; 1]);
    productX = product.addOrigin('x', @(x) x);
    productX2 = product.addOrigin('x^2', @(x) x(1,:).*x(2,:));

    n.addNode(product);
    n.addConnection(integratorX, productA);
    n.addConnection(integratorX, productB);
    twoDimProbe = n.addProbe(productX, 'output');
    productProbe = n.addProbe(productX2, 'output');
    
    sg3 = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    sum = CosinePopulation([1], sg3, 'sum');
    sumA = sum.addTermination('input', tau, [1 .5]);
    sumX = sum.addOrigin('x', @(x) x);  

    n.addNode(sum);
    n.addConnection(productX, sumA);
    sumProbe = n.addProbe(sumX, 'output');
    
    tn = TestNode('testnode');
    sg3 = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    testoutput = CosinePopulation([1], sg3, 'testoutput');
    testoutputA = testoutput.addTermination('input', tau, 1);
    testoutputX = testoutput.addOrigin('x', @(x) x);
    n.addNode(tn);
    n.addNode(testoutput);    
    n.addConnection(productX2, tn.terminations{1});
    n.addConnection(tn.origins{1}, testoutputA);
    testoutputProbe = n.addProbe(testoutputX, 'output'); 
    
%     % this is too slow to tolerate 
%     subTestNode = TestNode('subtestnode');
%     subnetA = Network(.001); 
%     subnetA.name = 'subnetA';
%     subnetA.addNode(subTestNode);
%     subnetA.exposeOrigin(subTestNode.origins{1});
%     subnetA.exposeTermination(subTestNode.terminations{1});
%     n.addNode(subnetA);
%     n.addConnection(input.getOrigin(), subnetA.terminations{1});  
    
%     % this is also quite slow 
%     sg4 = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
%     subPop = CosinePopulation([1], sg4, 'subpop');
%     subPopA = subPop.addTermination('input', tau, 1);
%     subPopX = subPop.addOrigin('x', @(x) x);
%     subnetB = Network(.001); 
%     subnetB.name = 'subnetB';
%     subnetB.addNode(subPop);
%     subnetB.exposeOrigin(subPopX);
%     subnetB.exposeTermination(subPopA);    
%     n.addNode(subnetB); 
%     n.addConnection(input.getOrigin(), subnetB.terminations{1});    
%     subnetBProbe = n.addProbe(subnetB.origins{1}, 'output'); 
    
    n.setSimulationMode(ModeConfigurable.DIRECT_MODE);
    T = 5;
    
    n.run(0, T);    
    assert(isempty(spikeProbe.getHistory))
    [t, h] = sumProbe.getHistory();
    assert(abs(h(end) - 7.5) < .5);
    [t, h] = inputProbe.getHistory();
    assert(abs(h(end) - 1) < 1e-3) 
    [t, h] = outputProbe.getHistory();
    assert(abs(h(end) - T) < 1e-2)    

    [t, h] = testoutputProbe.getHistory();
    assert(abs(h(end) - 46) < 1) 
    
    [t, h] = twoDimProbe.getHistory();    
    assert(size(h, 1) == 2, 'Dimension of history should be 2')

    [t, h] = productProbe.getHistory();
    assert(abs(h(end) - 24) < 1) 
    
%     [t, h] = subnetBProbe.getHistory();
%     assert(abs(h(end) - 1) < 1e-2);
end

