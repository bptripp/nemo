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
end
