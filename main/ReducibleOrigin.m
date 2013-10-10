% The Origin of a connection that can be simulated in reduced form. 

% TODO: how to handle changing state dimension?

classdef ReducibleOrigin < Origin
    
    properties (SetAccess = private) 
        group = [];
    end
    
    properties (SetAccess = protected)
        decoders = [];
    end
    
    properties (Access = public) 
        reducedDim = [] %empty means not reduced 
    end
    
    methods (Access = public) 
        
        function ro = ReducedOrigin(name, group)
            outputDim = group.spikeGenerator.n;
            ro = ro@Origin(name, outputDim);
            
        end
        
        function setDecoders(ro, decoders) 
            decod
        end
        
        % This should normally be called by the associated Node before 
        % completion of run(...), in order to make Node output available 
        % to other Nodes. When running in DIRECT_MODE, there is no activity
        % and setX(...) should be called instead. 
        % 
        % time: Simulation time (at end of run step)
        % activity: Firing rates of neurons in the population, or spikes with integral 1 (o.dim by 1)
        function setActivity(ro, time, activity)
            if isempty(ro.reducedDim) 
                setOutput(do, time, activity);                
            else
                assert(~isempty(ro.decoders), 'setDecoders(...) must be called before this ReducedOrigin is used');
                assert(size(activity, 2) == 1, 'Activity should be a column vector');
                assert(size(activity, 1) == size(do.decoders, 2), 'Activity should be a vector of length %i, not %i', ...
                    size(do.decoders, 2), size(activity, 1));

                output = ro.decoders(1:ro.reducedDim,:) * activity;
                setOutput(do, time, output);                
            end
        end
        
    end
end