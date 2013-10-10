% A leaky-integrate-and-fire population spiking model. 
classdef LIFSpikeGenerator < SpikeGenerator

    properties (Constant)
        STATE_CURRENT = 'current';        
        STATE_POTENTIAL = 'potential';   
    end
    
    properties (Access = public)
        dt = []; 
        R = 1; 
        tauRef = []; 
        tauRC = []; 
        tauAdapt = [];
        Ginc = [];
        scales = []; 
        biases = []; 
        V0 = []; 
        Gad0 = [];
        noiseSD = [];
    end
    
    properties (SetAccess = private) 
        V = []; 
        Gad = [];
        lastSpike = []; 

        lastTime = [];
        lastCurrentState = [];        
        lastPotentialState = [];        
    end
    
    methods (Access = public)

        % dt: internal time step (may be shorter than Network time step)
        % tauRef: refractory period (s) shared by all neurons
        % tauRC: membrane time constant (s) shared by all neurons
        % intercepts: TODO: how does this relate to radius?
        % maxRates: maximum firing rate of each neuron within represented
        %   region (1 x n)
        % V0: starting membrane potential of each neuron (1 x n or single
        %   value if shared)
        % tauAdapt (optional): 
        % Ginc (optional): 
        % Gad0 (optional): 
        % noiseSD (optional): standard deviation of intrinsic noise added
        %   each ms. If the spike generator is run in steps of a different 
        %   size, the noise amplitude is adjusted as if it is uncorrelated
        %   over time. 
        function sg = LIFSpikeGenerator(dt, tauRef, tauRC, intercepts, maxRates, V0, varargin)
            sg = sg@SpikeGenerator(length(intercepts));
            
            assert(dt < tauRef, 'Timestep must be smaller than spike refractory period');
            assert(length(intercepts) == length(maxRates), 'Expected one intercept and one max rate per neuron')
            assert(length(tauRef) == 1, 'Expected a single refractory period for the population');
            assert(length(tauRC) == 1, 'Expected a single RC time constant for the population');
            assert(length(V0) == length(intercepts) || length(V0) == 1, 'Expected either a single common V0 or one per neuron')
            assert(size(intercepts, 1) == 1, 'Expected row vectors for maxRates, intercepts, V0')
            assert(size(maxRates, 1) == 1, 'Expected row vectors for maxRates, intercepts, V0')
            assert(size(V0, 1) == 1, 'Expected row vectors for maxRates, intercepts, V0')

            sg.dt = dt;
            sg.tauRef = tauRef;
            sg.tauRC = tauRC;
            x = 1 ./ (1 - exp( (tauRef - (1 ./ maxRates)) / tauRC));
            sg.scales = (x - 1) ./ (1 - intercepts);
            sg.biases = 1 - sg.scales .* intercepts
            if length(V0) == 1
                sg.V0 = ones(size(sg.scales)) * V0;
            else 
                sg.V0 = V0;    
            end
            sg.tauAdapt = NemoUtils.getOptionalArg(varargin, 1, [], 'tauAdapt', size(sg.V0));
            sg.Ginc = NemoUtils.getOptionalArg(varargin, 2, zeros(size(sg.V0)), 'Ginc', size(sg.V0));
            sg.Gad0 = NemoUtils.getOptionalArg(varargin, 3, zeros(size(sg.V0)), 'Gad0', size(sg.V0));
            sg.noiseSD = NemoUtils.getOptionalArg(varargin, 4, [], 'noiseSD', 1);
            
            reset(sg);
        end
        
        % instances (optional): Number of instances to run in parallel
        %   (default 1). This is useful for vectorized simulation. If >1,
        %   should be called again with 1 before running within a normal
        %   simulation. 
        function reset(sg, varargin)
            reset@SpikeGenerator(sg);
            instances = NemoUtils.getOptionalArg(varargin, 1, 1, 'instances', 1);
            randomize = NemoUtils.getOptionalArg(varargin, 2, 0, 'randomize', 1);
            
            sg.V = repmat(sg.V0, instances, 1);
            if randomize
                sg.V = rand(size(sg.V));
            end
            sg.Gad = repmat(sg.Gad0, instances, 1);
            sg.lastSpike = -sg.tauRef*ones(size(sg.V));
        end
        
        % See SpikeGenerator.getRates(...)
        function rates = getRates(sg, drive, startTime, endTime)
            if endTime > startTime % run with noise
                T = endTime - startTime;
                dt = .001;
                spikes = zeros(size(drive));
                reset(sg, size(drive, 2), 1);
                for i = 1:floor(T/dt)
                    i
                    spikes = spikes + integrate(sg, drive, (i-1)*dt, i*dt);
                end
                spikes = spikes + integrate(sg, drive, floor(T/dt)*dt, endTime);
                rates = spikes / T;
            else 
                foo = ones(1, size(drive, 2));
                current = (sg.scales'*foo) .* drive + (sg.biases'*foo); 
                rates = zeros(size(current));
                indices = find(current > 1);
                rates(indices) = 1 ./ (sg.tauRef - sg.tauRC * log(1 - 1./current(indices)));    
            end
            reset(sg);
        end
        
        % See SpikeGenerator.integrate(...)
        function spikes = integrate(sg, drive, startTime, endTime)
            T = endTime - startTime;
            drive = addNoise(drive, sg.noiseSD, T);
            
            nSteps = ceil(T / sg.dt);
            stepLength = T / nSteps;
            spikes = zeros(length(sg.scales), size(sg.V, 1));
            for i = 1:nSteps
                spikes = spikes + step(sg, drive, startTime+(i-1)*stepLength, startTime+i*stepLength);
            end
            spikes = spikes > 0;
        end
        
        % see Probeable
        function names = getStateNames(sg)
            names = {SpikeGenerator.STATE_DRIVE, SpikeGenerator.STATE_SPIKES, SpikeGenerator.STATE_CURRENT, SpikeGenerator.STATE_POTENTIAL};
        end
        
        function dim = getDimension(sg, name)
            if strcmp(name, LIFSpikeGenerator.STATE_CURRENT) || strcmp(name, LIFSpikeGenerator.STATE_POTENTIAL)
                dim = sg.n;
            else 
                dim = getDimension@SpikeGenerator(sg, name);
            end            
        end        
        
        % see Probeable
        function [time, state] = getState(sg, name)
            if strcmp(name, LIFSpikeGenerator.STATE_CURRENT)
                time = sg.lastTime;
                state = sg.lastCurrentState;
            elseif strcmp(name, LIFSpikeGenerator.STATE_POTENTIAL)
                time = sg.lastTime;
                state = sg.lastPotentialState;
            else 
                [time, state] = getState@SpikeGenerator(sg, name);
            end
        end        
    end
    
    methods (Access = private) 
        function spikes = step(sg, drive, startTime, endTime)            
            assert(max(max(sg.lastSpike)) <= startTime, ...
                'Reset the SpikeGenerator before running from an earlier time (start time %f, last spike %f)', startTime, max(sg.lastSpike));
            instances = size(sg.V, 1);
            
            intTime = endTime - startTime;
            intTimes = min(intTime, max(0, endTime - sg.lastSpike - sg.tauRef));

            sg.lastTime = endTime;
            
            current = repmat(sg.biases, instances, 1) + repmat(sg.scales, instances, 1) .* drive';
            sg.lastCurrentState = current;
            
            dV = (1 / sg.tauRC) .* (sg.R*current - sg.V .* (1 + sg.R*sg.Gad)); 
            sg.V = sg.V + intTimes .* dV;
            sg.V = max(0, sg.V);
            sg.lastPotentialState = sg.V;
            
            if ~isempty(sg.tauAdapt)
                dGad = -sg.Gad ./ sg.tauAdapt; 
                sg.Gad = sg.Gad + intTimes .* dGad; 
            end
            
            spikeIndices = find(sg.V >= 1);
            sg.lastSpike(spikeIndices) = endTime - (sg.V(spikeIndices)-1) ./ dV(spikeIndices);
            sg.V(spikeIndices) = 0; 
            repGinc = repmat(sg.Ginc, instances, 1);
            sg.Gad(spikeIndices) = sg.Gad(spikeIndices) + repGinc(spikeIndices);  
            spikes = zeros(fliplr(size(drive)));
            spikes(spikeIndices) = 1; 
            spikes = spikes'; %TODO: API has drive and spikes transposed ... clean this up
        end
    end
end

function drive = addNoise(drive, noiseSD, T) 
    if ~isempty(noiseSD)
        SD = noiseSD / (T/.001)^.5;
        drive = drive + SD*randn(size(drive));
    end
end
