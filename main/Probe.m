% Reads state variables from Probeable objects (eg membrane potential from 
% a Neuron). Collected data can be displayed during a simluation or kept 
% for plotting afterwards.
% 
% Note that a Probe is typically added to a Network. The Network instructs
% the Probe to collect data every time step. 
classdef Probe < handle 

    properties (SetAccess = private)
        stateName = [];
    end
    
    properties (Access = private) 
        probeable = [];
        
        history = [];
        times = [];
        dim = [];
        len = [];
        
        spike = 0;
    end
    
    methods (Access = public)
        
        % probeable: The Probeable for which to record state history
        % stateName: Name of the state variable to collect
        function p = Probe(probeable, stateName)
            p.probeable = probeable;
            p.stateName = stateName;
            p.dim = probeable.getDimension(stateName);

            % check state to see if it consists of spikes ...
            [time, sample] = probeable.getState(stateName);
            if islogical(sample)
                p.spike = 1;
            end
            
            reset(p);
        end

        % Obtains and stores current state from the associated Probeable
        % 
        % t: time of collection (used as data time if Probeable returns NaN
        %   time)
        function collect(p, t)
            [time, state] = p.probeable.getState(p.stateName);
            if isnan(time)
                time = t;
            end
            assert(length(state) == p.dim, 'Expected state of dimension %i, not %i', p.dim, length(state));
            if ~isempty(p.times)
                lastPreviousTime = nanmax(p.times);
                assert(isnan(lastPreviousTime) || time > lastPreviousTime, 'Reset Probe before collecting data from earlier times. Has this Population been added to the Network?')
            end
            
            p.len = p.len + 1;
            if p.len > size(p.history, 2)
                increment = 1000;
                p.history = copy(p.history, increment, p.spike);
                p.times = copy(p.times, increment, 0);
            end
            
            if (p.spike)
                p.history(:,p.len) = logical(state);
            else 
                p.history(:,p.len) = state;
            end
            
            p.times(:,p.len) = time;
        end
        
        % tau (optional input): time constant of first-order filter to 
        %   apply to history (none if omitted)
        % time: time points for which history is available
        % history: collected history at corresponding time points
        function [time, history] = getHistory(p, varargin)
            %TODO: use NemoUtils.filter()
            time = p.times(:,1:p.len);
            history = p.history(:,1:p.len);
            
            if ~islogical(history) && ~isempty(varargin)
                tau = varargin{1};  
                [d, l] = size(history);
                history(:,1) = zeros(d, 1);
                for i = 2:l % filter
                    history(:,i) = history(:,i-1) + (history(:,i) - history(:,i-1)) * (time(i)-time(i-1)) / tau;
                end
            end
        end
        
        % times: Time stamps corresponding to each point in the state
        %   history
        function times = getTimes(p)
            times = p.times(:,1:p.len);
        end
        
        % Clears collected state history
        function reset(p)
            if (p.spike)
                p.history = false(p.dim, 1000);
            else 
                p.history = zeros(p.dim, 1000);
            end
            p.times = NaN*zeros(1, 1000);
            p.len = 0;
        end
    end
end

function big = copy(small, increment, spikes)
    if spikes
        big = false(size(small, 1), size(small, 2) + increment);
    else 
        big = zeros(size(small, 1), size(small, 2) + increment);
    end
    big(:,1:size(small,2)) = small;
end