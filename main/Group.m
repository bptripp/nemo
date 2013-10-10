% A group of point neurons with current-gated synapses that have similar 
% connections to other groups of neurons. 
classdef Group < Node
    
    properties (SetAccess = private)        
        radii = [];        
        offsets = []; % middle of the represented region in each dimension
        spikeGenerator = [];        
    end
    
    methods (Access = public)

        % spikeGenerator: the spike generation model for the population
        % name: Name of the population (should be unique within a Network)
        function p = Population(spikeGenerator, name)
            p.spikeGenerator = spikeGenerator;
            p.name = name;
        end

        % see Node
        function run(g, start, stop)
            drive = zeros(1, g.spikeGenerator.n);
            
            for i = 1:length(p.terminations)
                t = p.terminations{i};
                run(t, start, stop);
                drive = drive + getOutput(t);
            end

            simulationMode = p.simulationMode;
            spikeMode = (simulationMode == p.defaultMode);
            if simulationMode == p.constantRateMode
                start = stop;
            end

            activity = run(p.spikeGenerator, drive, start, stop, spikeMode); 
            for i = 1:length(p.origins)
                setActivity(p.origins{i}, stop, activity);
            end
        end
        
        function reset(p)
            reset@Node(p);
            reset(p.spikeGenerator);
        end

        % name: The name of the Origin
        % f: A function of the represented domain that is to be decoded by
        %   the new Origin
        % origin | T (optional): If a specific implementation is 
        %   needed, a DecodedOrigin object can be provided here. The caller 
        %   is responsible for calling the origin's findDecoders(..) method.
        %   Alternatively, this arg can be a positive scalar that indicates
        %   the run time to use for decoding (in seconds). The intrinsic
        %   noise in the spike generator is then accounted for in decoding.
        %   If T is not specified, by default intrisic noise is not used
        %   and external noise is added as a fraction of the maximum firing
        %   rate over the population.
        function origin = addOrigin(p, name, f, varargin)
            assert(strcmp(class(f), 'function_handle'), 'Expected a function handle');
            
            if isempty(varargin) || isempty(varargin{1})
                origin = DecodedOrigin(name, f, p);
                origin.findDecoders([], .002);
            elseif isa(varargin{1}, 'DecodedOrigin')
                origin = varargin{1};
            else 
                T = varargin{1};
                origin = DecodedOrigin(name, f, p);
                origin.findDecoders([], .002, T, 0);
            end
            
            p.origins{length(p.origins)+1} = origin;
        end

        % Adds an Origin that provides rate / spike outputs (depending on
        % the SimulationMode). The dimension of the Origin is the number of
        % neurons in the Population. 
        function origin = addAxonOrigin(p)
            n = p.spikeGenerator.n;
            origin = DecodedOrigin('AXON', @(x) zeros(n,size(x,2)), p);
            setDecoders(origin, eye(n));
            p.origins{length(p.origins)+1} = origin;
        end

        % Removes named Origin
        % name: Name of Origin to remove
        function removeOrigin(p, name)
            toRemove = [];
            for i = 1:length(p.origins)
                if strcmp(name, p.origins{i}.name)
                    toRemove = [toRemove i];
                end
            end
            toKeep = setdiff(1:length(p.origins), toRemove);
            p.origins = p.origins(toKeep);
        end
        
        % Removes all origins
        function removeOrigins(p)
            p.origins = {};
        end

        % name: Name of Termination (unique for the Population)
        % timeConstant: Time constant of a first-order filter through which
        %   inputs pass
        % transform: matrix that defines map from the input to the space 
        %   represented by the associated Population
        function termination = addTermination(p, name, timeConstant, transform)
            assert(size(transform, 1) == p.dimension, 'Transform should have %i rows', p.dimension);
            termination = Termination(name, timeConstant, transform, 1); %TODO: how to handle step ratios?
            p.terminations{length(p.terminations)+1} = termination;
        end
        
        % Removes named Termination
        % name: Name of Termination to remove
        function removeTermination(p, name)
            toRemove = [];
            for i = 1:length(p.terminations)
                if strcmp(name, p.terminations{i}.name)
                    toRemove = [toRemove i];
                end
            end
            toKeep = setdiff(1:length(p.terminations), toRemove);
            p.terminations = p.terminations(toKeep);
        end

        % Removes all Terminations
        function removeTerminations(p)
            p.terminations = {};
        end

    end    
    
end