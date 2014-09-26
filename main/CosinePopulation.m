% A population of neurons with cosine tuning.  
classdef CosinePopulation < Population
    properties (Access = public)
        encoders = [];
    end
    
    methods (Access = public) 
                
        function cp = CosinePopulation(radii, spikeGenerator, name, varargin)
            % radii: max amplitude of represented vector in each dimension
            %   beyond which saturation occurs, or radii of the ellipse 
            %   in the represented space that is accurately represented 
            % spikeGenerator: the spike generation model for the population
            % name: Name of the population (should be unique within a Network)
            % encoders (optional): length(radii) X n matrix of encoding
            %   vectors (default random uniformly distributed)
            % ellipsoidRegion (optional): 1 (default) if the represented region
            %   region is ellipsoidal; 0 if it is a box. 
            % offsets (optional): centre of point distribution in each
            %   dimension 
            
            cp = cp@Population(radii, spikeGenerator, name, varargin{2:end});
            n = spikeGenerator.n;
            
            if ~isempty(varargin) && ~isempty(varargin{1})
                cp.encoders = varargin{1};
                assert(size(cp.encoders, 2) == n, 'Expected same number of encoders as neurons in spike generator')
                assert(size(cp.encoders, 1) == length(radii), 'Expected encoders to have same length as radii')
            else 
                cp.encoders = Population.genRandomPoints(n, 1./radii, 1);
            end
        end
        
        % see Population.getDrive(...)
        function drive = getDrive(cp, x)
            drive = cp.encoders' * x;
        end
        
    end
    
    methods (Static)
        
        % Creates an n-dimensional population by cloning a 1D population. 
        % 
        % pop: population to replicate including origins (not terminations)
        function result = makeClonedPopulation(pop, n)
            assert(length(pop.radii) == 1, 'only 1D populations supported')

            nn = pop.spikeGenerator.n; % number of neurons in pop
            sg = CosinePopulation.makeAggregateSpikeGenerator(pop.spikeGenerator, n);   
            allEncoders = zeros(n, n*nn);
            for i = 1:n
                allEncoders(i,(i-1)*nn+(1:nn)) = pop.encoders;
            end
            result = CosinePopulation(repmat(pop.radii, n, 1), sg, pop.name, allEncoders);
            
            % add origins
            for i = 1:length(pop.origins)
                assert(isa(pop.origins{i}, 'DecodedOrigin'), 'only DecodedOrigins are supported');
                fun = @(x) zeros(n,size(x,2));
                do = DecodedOrigin(pop.origins{i}.name, fun, result); 
                decoders = pop.origins{i}.decoders;
                allDecoders = zeros(n, n*nn);
                for j = 1:n
                    allDecoders(j, (j-1)*nn+(1:nn)) = decoders;
                end
                setDecoders(do, allDecoders);
                result.addOrigin(pop.origins{i}.name, fun, do);
            end
        end
        
        % Note that only non-adapting LIF neurons are implemented. 
        function result = makeAggregateSpikeGenerator(sg, n)
            assert(isa(sg, 'LIFSpikeGenerator'), 'only LIF spike generators supported')            
            nn = length(sg.scales);
            result = LIFSpikeGenerator(sg.dt, sg.tauRef, sg.tauRC, zeros(1, nn*n), ones(1, nn*n), repmat(sg.V0, 1, n));
            result.scales = repmat(sg.scales, 1, n);
            result.biases = repmat(sg.biases, 1, n);
        end
        
        % Uses default encoders. 
        % 
        % reps: repetitions of low-D structure to aggregate
        % n: number of neurons per repetition
        % radii: low-D radii
        % tauRef: LIF refractory time
        % tauRC: LIF membrane time constant
        % intercepts: For all neurons (1 by reps*n) 
        % maxRates: For all neurons (1 by reps*n)
        % originFunctions: a struct with names and functions
        % 
        % Example: 
        % dim = 50; test = CosinePopulation.makeAggregatePopulation('product', dim, 100, 1, .002, .02, -1+2*rand(1,100*dim), 50+50*rand(1,100*dim), struct('x', @(x) x));
        function result = makeAggregatePopulation(name, reps, n, radii, tauRef, tauRC, intercepts, maxRates, originFunctions, varargin)
            assert(size(intercepts, 1) == 1, 'Intercepts should be 1 by reps*n');
            assert(size(intercepts, 2) == reps*n, 'Intercepts should be 1 by reps*n');
            
            encoders = [];
            if ~isempty(varargin)
                encoders = varargin{1};
                assert(size(encoders, 1) == length(radii), 'Encoders should be length(radii) by reps*n');
                assert(size(encoders, 2) == reps*n, 'Encoders should be length(radii) by reps*n');                
            end
            
            dt = .001;
            originNames = fieldnames(originFunctions);
            
            % make individual subpopulations to collapse later
            subpops = cell(reps,1);
            for i = 1:reps
                ind = (1:n)+(i-1)*n;
                sg = LIFSpikeGenerator(dt, tauRef, tauRC, intercepts(ind), maxRates(ind), 0);
                
                if isempty(encoders)
                    pop = CosinePopulation(radii, sg, name);
                else
                    pop = CosinePopulation(radii, sg, name, encoders(:,ind));
                end
                
                for j = 1:length(originNames)
                    pop.addOrigin(originNames{j}, originFunctions.(originNames{j}));
                end
                subpops{i} = pop;
            end
                        
            sg = LIFSpikeGenerator(dt, tauRef, tauRC, intercepts, maxRates, 0);
            allEncoders = zeros(length(radii)*reps, reps*n);
            d = length(radii);
            for i = 1:reps
                allEncoders((i-1)*d+(1:d), (i-1)*n+(1:n)) = subpops{i}.encoders;
            end
            result = CosinePopulation(repmat(radii, reps, 1), sg, name, allEncoders);
            
            for i = 1:length(originNames)
                of = originFunctions.(originNames{i});
                od = length(of(radii));
                fun = @(x) zeros(od*reps,size(x,2));  
                do = DecodedOrigin(originNames{i}, fun, result); 

                allDecoders = zeros(od*reps, reps*n);
                for j = 1:reps
                    allDecoders((j-1)*od+(1:od), (j-1)*n+(1:n)) = subpops{j}.getOrigin(originNames{i}).decoders;
                end
                setDecoders(do, allDecoders);
                result.addOrigin(originNames{i}, fun, do);
            end     
        end
        
    end
end

