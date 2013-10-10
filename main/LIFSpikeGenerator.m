% A leaky-integrate-and-fire population spiking model. 
% 
% TODO: update probe for grouped neurons
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
            [sg.scales sg.biases] = getScaleBias(tauRef, tauRC, maxRates, intercepts);
%             x = 1 ./ (1 - exp( (tauRef - (1 ./ maxRates)) / tauRC));
%             sg.scales = (x - 1) ./ (1 - intercepts);
%             sg.biases = 1 - sg.scales .* intercepts;
            if length(V0) == 1
                sg.V0 = ones(size(sg.scales)) * V0;
            else 
                sg.V0 = V0;    
            end
            sg.tauAdapt = NemoUtils.getOptionalArg(varargin, 1, [], 'tauAdapt', size(sg.V0));
            sg.Ginc = NemoUtils.getOptionalArg(varargin, 2, zeros(size(sg.V0)), 'Ginc', size(sg.V0));
            sg.Gad0 = NemoUtils.getOptionalArg(varargin, 3, zeros(size(sg.V0)), 'Gad0', size(sg.V0));
            sg.noiseSD = NemoUtils.getOptionalArg(varargin, 4, [], 'noiseSD', 1);
            
            sg.mergedCount = ones(sg.n, 1);
            
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
        function rates = getRates(sg, drive, startTime, endTime, varargin)
            indices = NemoUtils.getOptionalArg(varargin, 1, 1:sg.n, 'indices', []);
            assert(length(indices) == size(drive, 1), 'indices and drive should be the same size');
            
            if endTime > startTime % run with noise
                T = endTime - startTime;
                dt = .001;
                spikes = zeros(size(drive));
                reset(sg, size(drive, 2), 1);
                for i = 1:floor(T/dt)
                    i
                    spikes = spikes + integrate(sg, drive, (i-1)*dt, i*dt, indices);
                end
                spikes = spikes + integrate(sg, drive, floor(T/dt)*dt, endTime, indices);
                rates = spikes / T;
            else 
                foo = ones(1, size(drive, 2));
                current = (sg.scales(indices)'*foo) .* drive + (sg.biases(indices)'*foo); 
                rates = zeros(size(current));
                activeIndices = find(current > 1);
                rates(activeIndices) = 1 ./ (sg.tauRef - sg.tauRC * log(1 - 1./current(activeIndices)));    
            end
            reset(sg);
        end        

        % See SpikeGenerator.addMerge(...)
        function addMerge(sg, indices)
            assert(min(indices) > 0, 'Indices must be >0');
            assert(max(indices) <= length(sg.mergedCount), sprintf('Indices must be <= %i', length(sg.mergedCount)));
            
            sg.mergedCount = [sg.mergedCount; sum(sg.mergedCount(indices))];
            
            totalCount = sum(sg.mergedCount(indices));
            
            [maxRates intercepts] = getMaxRateIntercept(sg.tauRef, sg.tauRC, sg.scales(indices), sg.biases(indices));
            maxRate = maxRates * sg.mergedCount(indices) / totalCount;
            intercept = intercepts * sg.mergedCount(indices) / totalCount;
            [scale bias] = getScaleBias(sg.tauRef, sg.tauRC, maxRate, intercept);
            sg.scales = [sg.scales scale];
            sg.biases = [sg.biases bias];            
            
%             sg.scales = [sg.scales sg.scales(indices) * sg.mergedCount(indices) / totalCount];
%             sg.biases = [sg.biases sg.biases(indices) * sg.mergedCount(indices) / totalCount];
            sg.V0 = [sg.V0 sg.V0(indices) * sg.mergedCount(indices) / totalCount];

            sg.Ginc = [sg.Ginc sg.Ginc(indices) * sg.mergedCount(indices) / totalCount];
            sg.Gad0 = [sg.Gad0 sg.Gad0(indices) * sg.mergedCount(indices) / totalCount]; 
            if ~isempty(sg.tauAdapt)
                sg.tauAdapt = [sg.tauAdapt sg.tauAdapt(indices) * sg.mergedCount(indices) / totalCount];
            end
            
            sg.lastSpike = [sg.lastSpike -sg.tauRef];
        end
        
        % See SpikeGenerator.removeMerge(...)
        function removeMerge(sg, indices)
            assert(min(indices) > sg.n, 'Indices must be >n (cannot removed non-merged neurons)');
            assert(max(indices) <= length(sg.mergedCount), sprintf('Indices must be <= %f', length(sg.mergedCount)));
            
            toKeep = setdiff(1:length(sg.mergedCount), indices);
            
            sg.mergedCount = sg.mergedCount(toKeep);
            sg.scales = sg.scales(toKeep);
            sg.biases = sg.biases(toKeep);
            sg.V0 = sg.V0(toKeep);
            
            sg.Ginc = sg.Ginc(toKeep);
            sg.Gad0 = sg.Gad0(toKeep);
            if ~isempty(sg.tauAdapt)
                sg.tauAdapt = sg.tauAdapt(toKeep);
            end
            
            sg.lastSpike = sg.lastSpike(toKeep);
        end
        
        % See SpikeGenerator.integrate(...)
        function spikes = integrate(sg, drive, startTime, endTime, varargin)
            indices = NemoUtils.getOptionalArg(varargin, 1, 1:sg.n, 'indices', []); 
            
            T = endTime - startTime;
            drive = addNoise(drive, sg.noiseSD, T);
            
            nSteps = ceil(T / sg.dt);
            stepLength = T / nSteps;
            spikes = zeros(length(indices), size(sg.V, 1));
            for i = 1:nSteps
                spikes = spikes + step(sg, drive, startTime+(i-1)*stepLength, startTime+i*stepLength, indices);
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
        function spikes = step(sg, drive, startTime, endTime, indices)  
            assert(max(max(sg.lastSpike)) <= startTime, ...
                'Reset the SpikeGenerator before running from an earlier time (start time %f, last spike %f)', startTime, max(sg.lastSpike));
            instances = size(sg.V, 1);
            
            intTime = endTime - startTime;
            intTimes = min(intTime, max(0, endTime - sg.lastSpike - sg.tauRef));

            sg.lastTime = endTime;
            
            current = repmat(sg.biases(indices), instances, 1) + repmat(sg.scales(indices), instances, 1) .* drive';
            sg.lastCurrentState(indices) = current;
            
            dV = (1 / sg.tauRC) .* (sg.R*current - sg.V(indices) .* (1 + sg.R*sg.Gad(indices))); 
            sg.V(indices) = sg.V(indices) + intTimes(indices) .* dV;
            sg.V(indices) = max(0, sg.V(indices));
            sg.lastPotentialState(indices) = sg.V(indices);
            
            if ~isempty(sg.tauAdapt)
                dGad = -sg.Gad(indices) ./ sg.tauAdapt(indices); 
                sg.Gad(indices) = sg.Gad(indices) + intTimes .* dGad; 
            end
            
            % we can deal in spike indices below because there will
            % be no spikes outside indices ... 
            spikeIndices = find(sg.V >= 1);
                sg.lastSpike(spikeIndices) = endTime - (sg.V(spikeIndices)-1) ./ dV(ismember(indices, spikeIndices));
            sg.V(spikeIndices) = 0; 
            repGinc = repmat(sg.Ginc, instances, 1);
            sg.Gad(spikeIndices) = sg.Gad(spikeIndices) + repGinc(spikeIndices);  
            spikes = zeros(fliplr(size(sg.scales)));
            spikes(spikeIndices) = 1; 
%             spikes = spikes'; %TODO: API has drive and spikes transposed ... clean this up
            
            spikes = spikes(indices); % return only spikes from requested neurons
        end
    end
end

function drive = addNoise(drive, noiseSD, T) 
    if ~isempty(noiseSD)
        SD = noiseSD / (T/.001)^.5;
        drive = drive + SD*randn(size(drive));
    end
end

function [scales biases] = getScaleBias(tauRef, tauRC, maxRates, intercepts) 
    x = 1 ./ (1 - exp( (tauRef - (1 ./ maxRates)) / tauRC));
    scales = (x - 1) ./ (1 - intercepts);
    biases = 1 - scales .* intercepts;
end

function [maxRates intercepts] = getMaxRateIntercept(tauRef, tauRC, scales, biases)
    intercepts = (1-biases) ./ scales;
    maxRates = ( tauRef - tauRC * log(1 - 1./(scales + biases)) ).^-1;
end