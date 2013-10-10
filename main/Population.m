% A group of neurons that represent a variable together. 
classdef Population < Node
    
    properties (SetAccess = private)        
        radii = [];        
        offsets = []; % middle of the represented region in each dimension
        spikeGenerator = [];        
        dimension = [];
        ellipsoidRegion = [];
        indices = [];
    end
    
    properties (SetAccess = protected)
        Z = [];
    end
    
    properties (Access = private)
        % These are copied here because calling them directly every step
        % takes about twice as long (about 5% of run time in testInput)
        directMode = ModeConfigurable.DIRECT_MODE;
        defaultMode = ModeConfigurable.DEFAULT_MODE;
        constantRateMode = ModeConfigurable.CONSTANT_RATE_MODE;
    end

    methods (Abstract)
        % x: instantaneous variable represented by the population
        % indices (optional): indices of neurons and/or clusters for which
        %   to return drive (default 1:n)
        % drive: radial drive of each neuron (proportional to net driving
        %   current)
        drive = getDrive(p, x, varargin);        
    end
    
    methods (Access = public)

        % radii: max amplitude of represented vector in each dimension
        %   beyond which saturation occurs, or radii of the ellipse 
        %   in the represented space that is accurately represented 
        % spikeGenerator: the spike generation model for the population
        % name: Name of the population (should be unique within a Network)
        % ellipsoidRegion (optional): 1 (default) if the represented region
        %   region is ellipsoidal; 0 if it is a box. 
        % offsets (optional): centre of point distribution in each
        %   dimension 
        function p = Population(radii, spikeGenerator, name, varargin)
            assert(size(radii, 2) == 1, 'Expected column vector for radii')
            p.radii = radii;  
            p.offsets = NemoUtils.getOptionalArg(varargin, 2, zeros(size(radii)), 'offsets', size(radii));
            p.spikeGenerator = spikeGenerator;
            p.dimension = length(radii);
            p.name = name;
            p.ellipsoidRegion = 1;
            if ~isempty(varargin)
                p.ellipsoidRegion = varargin{1};
            end
        end
        
        % Finds a pattern of neuron clustering that minimizes simulation 
        % time within the error tolerances of this Population's
        % DecodedOrigins. Updates p.indices.
        function updateClusterIndices(p)
            n = p.spikeGenerator.n;
            
            %TODO: upper and lower limits 
            %TODO: interpolate from selected levels
            %TODO: bidirectional search from arbitrary point (importantly,
            %   the last point in a simulation where salience changes
            %   gradually -- use optional arg indices)
            
            doInd = []; % indices of decoded origins within p.origins            
            outInd = []; % start and end indices of each decoded origin's outputs within error vector
            for i = 1:length(p.origins)
                if isa(p.origins{i}, 'DecodedOrigin')
                    doInd = [doInd; i];   
                    
                    highestExistingInd = 0;
                    if ~isempty(outInd) 
                        highestExistingInd = max(outInd(:,2));
                    end
                    
                    outInd = [outInd; highestExistingInd+1 highestExistingInd+p.origins{i}.dim]; 
                end
            end
            
            tolerances = zeros(1, max(outInd(:,2)));
            for i = 1:length(doInd)
                tolerances(outInd(i,:)) = p.origins{doInd(i)}.tolerance;
            end
            
            clusterInd = 2*n-1; %start at root cluster 
            currentErrors = concatErrors(clusterInd, outInd, p.origins(doInd));
            while max(currentErrors > tolerances)
                values = getSplitValues(clusterInd, n, outInd, p.origins(doInd));
                correspondance = values * max(0, (currentErrors - tolerances)');
                [Y, I] = max(correspondance);
                clusterInd = [clusterInd(1:I-1) clusterInd(I+1:end) p.Z(clusterInd(I)-n,:)];
                currentErrors = concatErrors(clusterInd, outInd, p.origins(doInd));
            end
            
            p.indices = clusterInd;
        end
        
        % see Node
        % TODO: use current cluster indices in correct mode
        function run(p, start, stop)
            x = zeros(size(p.radii));
            biasDrive = zeros(1, p.spikeGenerator.n);
            for i = 1:length(p.terminations)
                %TODO: run with weights
                %TODO: transition between weights and encoders over time
                t = p.terminations{i};
                run(t, start, stop);
                if ~isempty(t.biasEncoders)
                    biasDrive = biasDrive + getOutput(t) * t.biasEncoders;
                else 
                    x = x + getOutput(t);
                end
            end

            simulationMode = p.simulationMode;
            if (simulationMode == p.directMode) 
                for i = 1:length(p.origins)
                    setX(p.origins{i}, stop, x);
                end
            else  
                drive = getDrive(p, x) + biasDrive';
            
                spikeMode = (simulationMode == p.defaultMode);
                if simulationMode == p.constantRateMode
                    start = stop;
                end
                
                activity = run(p.spikeGenerator, drive, start, stop, spikeMode); 
                for i = 1:length(p.origins)
                    setActivity(p.origins{i}, stop, activity);
                end
            end
        end
        
        function reset(p)
            reset@Node(p);
            reset(p.spikeGenerator);
        end

        % domain: limits of represented area in each dimension
        function domain = getDomain(p)
            domain = p.offsets + [-p.radii p.radii];
        end
        
        % dim: dimension of the represented variable
        function dim = getDimension(p)
            dim = length(p.radii);
        end
        
        % x: value of represented vector
        % startTime: see SpikeGenerator.getRates(...)
        % endTime: see SpikeGenerator.getRates(...)
        % indices (optional): indices of neurons and/or clusters for which
        %   to return rates (default 1:n)
        % rates: spike rates of neurons in the population (n by 1) at endTime
        function rates = getRates(p, x, startTime, endTime, varargin)
            indices = NemoUtils.getOptionalArg(varargin, 1, 1:p.spikeGenerator.n, 'indices', []);
            drive = getDrive(p, x, indices);
            rates = getRates(p.spikeGenerator, drive, startTime, endTime, indices);
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
                origin.findDecoders([], 0, 0, .05);
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

        % Adds a termination for bias signals associated with a Parisien
        % transform. 
        % 
        % name: Name of Termination (unique for the Population)
        % timeConstant: Time constant of a first-order filter through which
        %   inputs pass
        % transform: matrix that defines map from the input to the space 
        %   represented by the associated Population
        function termination = addBiasTermination(p, name, timeConstant, transform, bias, biasEncoders)
            termination = Termination(name, timeConstant, transform, 1);
            if ~isempty(bias)
                setBias(termination, bias);
            end
            setBiasEncoders(termination, biasEncoders);
            p.terminations{length(p.terminations)+1} = termination;
        end
        
    end    
    
    methods (Static) 
        % n: number of points to generate
        % radii: radii of ellipsoidal region in which points should lie
        %   (column vector)
        % surface: 0 if points should lie within region with uniform
        %   density, 1 if they should lie on the surface, 2 if sampled from
        %   a Gaussian (only valid with ellipsoid, below)
        % ellipsoid (optional): 1 (default) if points should be sampled in 
        %   an ellipsoid; 0 if in a box 
        % offsets (optional): centre of point distribution in each
        %   dimension
        % points: randomly-selected points (length(radii) x n)
        function points = genRandomPoints(n, radii, surface, varargin)
            offsets = NemoUtils.getOptionalArg(varargin, 2, zeros(size(radii)), 'offsets', size(radii));
            assert(size(radii, 2) == 1, 'Radii should be a column vector')
            dim = length(radii);

            box = 0;
            if ~isempty(varargin)
                box = ~varargin{1};
            end
            
            if box
                u = rand(length(radii), n);
                if surface
                    % push each point to edge along a random dimension ... 
                    indices = sub2ind([dim n], randi(dim, 1, n), 1:n);
                    u(indices) = u(indices) > .5;
                end
                points = -radii*ones(1,n) + 2*radii*ones(1,n) .* u;
            else 
                directions = randn(dim, n);
                                
                if surface == 2
                    points = diag(radii) * directions;
                else                     
                    len = sum(directions.^2, 1).^(1/2);
                    points = directions ./ (1./radii * len);
                    if surface == 0
                        points = points .* (ones(size(radii)) * rand(1,n).^(1/dim));    
                    end    
                end
            end
            points = offsets * ones(1, size(points,2)) + points;
        end

    end
end

% Concatenates errors from DecodedOrigins
function errors = concatErrors(indices, outInd, origins)
    errors = zeros(1, max(outInd(:,2)));
    for i = 1:length(origins)
        errors(outInd(i,:)) = getErrorEstimate(origins{i}, indices);
    end
end

function values = getSplitValues(clusterInd, n, outInd, origins)
    splittables = find(clusterInd > n);
    values = zeros(length(clusterInd), max(outInd(:,2)));
    for i = 1:length(origins)
        values(splittables,outInd(i,:)) = origins{i}.splitValues(:,clusterInd(splittables)-n)';
    end
end

