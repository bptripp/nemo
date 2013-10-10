% Unit tests for PopulationModeModel

function testPopulationMode() 
    showPlots = 1;

    % neuron parameters ... 
    tauRef = .0025;
    tauRC = .02;
    n = 200;
    intercepts = -1+2*rand(1,n);
    maxRates = 80+50*rand(1,n);
    V0 = 0;
    sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    
%     % create 1D population ... 
%     pop1 = CosinePopulation([1], sg, 'test1');  
%     pop1.addOrigin('x', @(x) x);
%     pop1.addOrigin('x2', @(x) [x.^2; abs(x)]);
%     initPopulationMode(pop1);
%     testNoiseGeneration(pop1);
%     test1DSimulation(pop1, showPlots)
%     
%     % create 2D population ... 
%     pop2 = CosinePopulation([1; 1], sg, 'test2');  
%     pop2.addOrigin('x', @(x) x);
%     initPopulationMode(pop2);
%     testSimulation(pop2, showPlots)
%     
%     % create 3D population ... 
%     pop3 = CosinePopulation([1; 1; 1], sg, 'test3'); 
%     pop3.addOrigin('x', @(x) x);    
%     initPopulationMode(pop3);
%     testSimulation(pop3, showPlots)
    
    % create 25D population ... 
    pop25 = CosinePopulation(ones(25, 1), sg, 'test25');  
    pop25.addOrigin('x', @(x) x);
    initPopulationMode(pop25); 
    testSimulation(pop25, showPlots)
end

function testNoiseGeneration(pop)
    modelNoise = generateNoise(pop.populationModeModel, 1000);
    simNoise = getNoiseSamples(pop.populationModeModel, .5, .001, 1); % 0.5 should have roughly average properties
    varRatios = var(modelNoise') ./ var(simNoise');
    assert(abs(mean(varRatios) - 1) < 1/4, sprintf('Noise variance ratios are not close to one: %s', sprintf('%f ', varRatios)))
    
    % TODO: test correlation, autocorrelation
end

function testSimulation(pop, showPlots)
    input = FunctionInput(@(t) t);

    tau = .01;
    popA = pop.addTermination('input', tau, ones(length(pop.radii), 1));

    n = Network(.001);
    n.addNode(input);
    n.addNode(pop);
    n.addConnection(input.origins{1}, popA);
    outputProbe = n.addProbe(getOrigin(pop, 'x'), 'output');

    T = 3;
    n.reset();
    setSimulationMode(pop, ModeConfigurable.DEFAULT_MODE);
    n.run(0, T);
    
    [time, defaultModeX] = getHistory(outputProbe);
    if showPlots
        plotProbe(outputProbe, .01)
    end
    
    n.reset();
    setSimulationMode(pop, ModeConfigurable.POPULATION_MODE);
    n.run(0, T);
    
    [time, populationModeX] = getHistory(outputProbe);
    if showPlots
        plotProbe(outputProbe, .01)
    end
    
    ind = length(time)+(-20:0);
    assert(abs(mean(defaultModeX(ind)) - mean(populationModeX(ind))) < .1) 
end

% This is mostly a duplicate of testSimulation, but with an additional 
% origin with nonlinear functions. 
function test1DSimulation(pop, showPlots)
    input = FunctionInput(@(t) t);

    tau = .01;
    popA = pop.addTermination('input', tau, 1);

    n = Network(.001);
    n.addNode(input);
    n.addNode(pop);
    n.addConnection(input.origins{1}, popA);
    outputProbe = n.addProbe(getOrigin(pop, 'x'), 'output');
    output2Probe = n.addProbe(getOrigin(pop, 'x2'), 'output');

    T = 3;
    n.reset();
    setSimulationMode(pop, ModeConfigurable.DEFAULT_MODE);
    n.run(0, T);
    
    [time, defaultModeX] = getHistory(outputProbe);
    [time, defaultModeX2] = getHistory(output2Probe);
    if showPlots
        plotProbe(outputProbe, .01)
        plotProbe(output2Probe, .01)
    end
    
    n.reset();
    setSimulationMode(pop, ModeConfigurable.POPULATION_MODE);
    n.run(0, T);
    
    [time, populationModeX] = getHistory(outputProbe);
    [time, populationModeX2] = getHistory(output2Probe);
    if showPlots
        plotProbe(outputProbe, .01)
        plotProbe(output2Probe, .01)
    end
    
    ind = length(time)+(-20:0);
    assert(abs(mean(defaultModeX(ind)) - mean(populationModeX(ind))) < .1) 
    assert(abs(mean(defaultModeX2(1,ind)) - mean(populationModeX2(1,ind))) < .1) 
    assert(abs(mean(defaultModeX2(2,ind)) - mean(populationModeX2(2,ind))) < .1) 
end
