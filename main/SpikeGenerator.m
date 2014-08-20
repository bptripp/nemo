% A spike generation model for a population. (Spiking for the whole
% population is simulated together to allow vectorization.)
classdef SpikeGenerator < Probeable 
    
    properties (Constant)
        STATE_DRIVE = 'drive';
        STATE_SPIKES = 'spikes';
    end
    
    properties (Access = private)
        lastTime = [];
        lastDriveState = [];
        lastSpikeState = [];
    end
    
    properties (SetAccess = private)
        n = 0;
    end
    
    methods (Abstract)

        % drive: net radial drive entering each neuron's spike
        %   generator (n by 1)
        % startTime: beginning of simulation time (relevant if rate is not
        %   a static function of current, eg with adaptation)
        % endTime: end of simulation time (if a static rate is needed, this
        %   may be the same as startTime)
        % rates: spike rates (n by 1) at endTime
        rates = getRates(sg, drive, startTime, endTime);
        
        % drive: net radial drive entering each neuron's spike
        %   generator (n by 1)
        % startTime: beginning of simulation time 
        % endTime: end of simulation time
        % spikes: 1 for each neuron that spiked during the integration
        %   step, 0 for those that didn't (n by 1)
        spikes = integrate(sg, drive, startTime, endTime);
    end
    
    methods (Access = public)

        % n: number of neurons in the population
        function sg = SpikeGenerator(n)
            sg.n = n;
            sg.reset();
        end

        % see Probeable
        function names = getStateNames(sg)
            names = {SpikeGenerator.STATE_DRIVE, SpikeGenerator.STATE_SPIKES};
        end

        % see Probeable
        function dim = getDimension(sg, name)
            knownState = strcmp(name, SpikeGenerator.STATE_DRIVE) || strcmp(name, SpikeGenerator.STATE_SPIKES);
            assert (knownState, 'State %s unknown', name);
            dim = sg.n;
        end

        % see Probeable
        function [time, state] = getState(sg, name)
            isDriveState = strcmp(name, SpikeGenerator.STATE_DRIVE);
            isSpikeState = strcmp(name, SpikeGenerator.STATE_SPIKES);
            assert (isDriveState || isSpikeState, 'State %s unknown', name);
            time = sg.lastTime;
            if isDriveState
                state = sg.lastDriveState;
            else 
                state = sg.lastSpikeState;
            end
        end
        
        function reset(sg)
            sg.lastTime = NaN;
            sg.lastDriveState = zeros(sg.n,1);
            sg.lastSpikeState = false(sg.n,1);
        end

        % drive: net radial drive entering each neuron's spike
        %   generator (n by 1)
        % startTime: beginning of simulation time 
        % endTime: end of simulation time
        % output: spike impulses (integral one) or rates depending on
        %   the spikeMode property (n by 1)
        function output = run(sg, drive, startTime, endTime, spikeMode)
            sg.lastTime = endTime;
            if spikeMode
                output = integrate(sg, drive, startTime, endTime) / (endTime-startTime);
                sg.lastDriveState = drive;
                sg.lastSpikeState = logical(output);
            else 
                output = getRates(sg, drive, startTime, startTime); %run at an instant, not for an interval
                sg.lastDriveState = drive;
                sg.lastSpikeState = false(size(output));
            end
        end
        
    end
end