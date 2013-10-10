% A network element, such as a neuron, population, or abstract input.  
classdef Node < ModeConfigurable
    
    properties (Access = public)        
        name = [];
    end
    
    properties (SetAccess = protected)
        origins = cell(1, 0);
        terminations = cell(1, 0);
    end
    
    methods (Abstract)
        
        % start: start time (s)
        % stop: stop time (s)
        % originNames: cell array of names of origins that must be updated
        %   in this run (empty = all)
        % originDims: cell array of vectors of dimension numbers of the
        %   corresponding originNames (empty = all dimensions)
        run(n, start, stop, originNames, originDims);
        
    end
    
    methods
        
        %TODO: reset to last time step? or another point? how much history
        %do we want to save? not much at anything other than coarse (could
        %be at regular intervals)

        % resets origins and terminations
        function reset(n)
            for i = 1:length(n.origins)
                reset(n.origins{i});
            end
            for i = 1:length(n.terminations)
                reset(n.terminations{i});
            end
        end
        
        % name: name of origin
        % origin: origin by the given name if it exists, otherwise []
        function origin = getOrigin(n, name)
            origin = [];
            for i = 1:length(n.origins)
                if strcmp(n.origins{i}.name, name)
                    origin = n.origins{i};
                    break;
                end
            end
        end
    end
end