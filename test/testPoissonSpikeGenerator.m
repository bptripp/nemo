%TODO: do I not need the correlation correction for spike counts with high
%rates?
function testPoissonSpikeGenerator()
    n = 10;
    intercepts = -1+2*rand(1,n);
    maxRates = 20+10*rand(1,n);
    centres = rand(1,n);
    f = @(x, y) max(.3-abs(x-y)/2,0);
    correlations = PoissonSpikeGenerator.genCorrelations(centres, f, .1);
%     correlations = eye(n);
    sg = PoissonSpikeGenerator(intercepts, maxRates, correlations);

    T = 1;
    dt = .001;
    time = 0:dt:T;
%     drive = -1:.4:1;
%     drives = ones(n,1) * drive;
%     rates = getRates(sg, drives, 0, T)
% 
%     integral = zeros(size(rates));
%     for i = 1:length(drive)  
%         sg.reset();
%         for j = 1:(length(time)-1)
%             output = sg.run(drives(:,i), time(j), time(j+1), 1);
%             integral(:,i) = integral(:,i) + output * (time(j+1)-time(j));
%         end        
%     end
%     
%     integral
%     figure, hold on
%     plot(drive, getRates(sg, drives, 0, 0))
%     plot(drive, rates', ':')
%     plot(drive, integral'/T, '--')
%     title('Ideal rates, count samples, and spike samples')
    
    drive = ones(n,1);
    
    spikes = zeros(n, length(time)-1);
    for i = 1:(length(time)-1)
        output = sg.run(drive, time(i), time(i+1), 1);
        spikes(:,i) = output;
    end
    figure, hold on
    for i = 1:n
        plot(find(spikes(i,:)), i, 'k.');
    end
    pause 
    
    figure, PoissonSpikeGenerator.plotCorrelationVsDistance(centres, correlations) 
    title('Specified Correlations')
    
    nTrials = 100;
    drive = ones(n,1);
    
    counts = zeros(n, nTrials);
    for i = 1:nTrials
        i
        counts(:,i) = getRates(sg, drive, 0, T);
    end
    rho = corr(counts');
    figure, PoissonSpikeGenerator.plotCorrelationVsDistance(centres, rho) 
    title('Count Correlations')

    spikeCounts = zeros(n, nTrials);
    for i = 1:nTrials
        i
        sg.reset();
        for j = 1:(length(time)-1)
            output = sg.run(drive, time(j), time(j+1), 1);
            spikeCounts(:,i) = spikeCounts(:,i) + output * (time(j+1)-time(j));
        end        
    end
    rho = corr(spikeCounts');
    figure, PoissonSpikeGenerator.plotCorrelationVsDistance(centres, rho) 
    title('Spike Count Correlations')
   
%     keyboard
end