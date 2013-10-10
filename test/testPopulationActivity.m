% Compares rates and activity in a population
% TODO: unit test with random input in multidimensional population
function testPopulationActivity() 
    tauRef = .0025;
    tauRC = .02;
    n = 10;
    intercepts = -1+2*rand(1,n);
    maxRates = 100+200*rand(1,n);
    V0 = 0;
    sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    test = CosinePopulation([1], sg, 'test');
    
    tau = .01;
    testInput = test.addTermination('input', tau, 1);
    ratesOrigin = test.addAxonOrigin()
    
    x = -1:.5:1;
    dt = .001;
    T = 1;
    time = 0:dt:T;
    rates = zeros(n, length(x));
    counts = zeros(n, length(x));
    for i = 1:length(x)
        i
        rates(:,i) = test.getRates(x(i), 0, 0);
        
        testInput.setInput(x(i));
        for j = 2:length(time)
            test.run(time(j-1), time(j));
            counts(:,i) = counts(:,i) + dt*ratesOrigin.getOutput();
        end
        reset(test);
    end
    figure, plot(x, rates'), title('population rates')
    figure, plot(x, rates', 'k', x, counts', 'r'), title('population spike counts per s (red) and rates (black)')
end
