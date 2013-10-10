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
        run(n, start, stop);
    end
    
    methods
        
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