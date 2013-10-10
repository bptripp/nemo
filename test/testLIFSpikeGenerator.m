% TODO: use unit testing framework
function testLIFSpikeGenerator() 
    tauRef = .0025;
    tauRC = .02;
    n = 50;
    intercepts = -1+2*rand(1,n);
    maxRates = 100+200*rand(1,n);
    V0 = 0;
    sg = LIFSpikeGenerator(.0005, tauRef, tauRC, intercepts, maxRates, V0);
    
    drive = -1:.4:1;
    drives = ones(n,1) * drive;
    rates = getRates(sg, drives, 0, 0);

    dt = .001;
    T = 1;
    time = 0:dt:T;
    integral = zeros(size(rates));
    for i = 1:length(drive)  
        sg.reset();
        for j = 1:(length(time)-1)
            output = sg.run(drives(:,i), time(j), time(j+1), 1);
            integral(:,i) = integral(:,i) + output * (time(j+1)-time(j));
        end        
    end
    
    integral
    figure, hold on
    plot(drive, rates')
    plot(drive, integral'/T, '--')
end