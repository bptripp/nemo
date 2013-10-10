% A population of neurons in which tuning is determined by user-defined
% functions (eg Gaussian; centre-surround)
classdef FunctionPopulation < Population 
    
    properties (Access = public)
        functions = [];
    end
    
    methods (Access = public)
        
        function fp = FunctionPopulation(radii, spikeGenerator, name, functions, varargin) 
            assert(spikeGenerator.n == length(functions), 'Expected same number of functions as size of spike generator');
            fp = fp@Population(radii, spikeGenerator, name, varargin{:});
            fp.functions = functions;
        end
        
        % see Population.getDrive(...)
        function drive = getDrive(fp, x)
            assert(size(x,1) == length(fp.radii), 'Expected column vectors');
            
            tic
            drive = zeros(length(fp.functions), size(x, 2));
            for i = 1:size(drive, 1)
                f = fp.functions{i};
                drive(i,:) = f(x);                
            end
            toc
        end
        
    end
        
end