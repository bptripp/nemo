% A population of neurons with Gaussian tuning (i.e. net drive is a 
% Gaussian function of the represented variable). 
classdef GaussianPopulation < Population
    
    properties (Access = public)
        centres = []; % dim by n
        covariances = []; % dim by dim by n
    end
    
    methods (Access = public)
        
        function gp = GaussianPopulation(radii, spikeGenerator, name, centres, covariances, varargin) 
            assert(size(centres, 1) == length(radii), 'Matrix of centres should have same number of rows as length of radii')
            assert(size(centres, 2) == spikeGenerator.n, 'Expected same number of centres as n of spikeGenerator')
            assert(size(covariances, 1) == size(centres, 1), 'Expected covariances with same dimension as centres')
            assert(size(covariances, 2) == size(centres, 1), 'Expected covariances with same dimension as centres')
            assert(size(covariances, 3) == size(centres, 2), 'Expected same number of covariances as centres')
            
            gp = gp@Population(radii, spikeGenerator, name, varargin{:});
            gp.centres = centres;
            gp.covariances = covariances;
        end
        
        % see Population.getDrive(...)
        function drive = getDrive(gp, x)
            [dimx, nx] = size(x);
            assert(dimx == length(gp.radii), 'Expected values with same dimension as the population');
            
            nn = size(gp.centres, 2); % # of neurons
            drive = zeros(nn, nx);
            for i = 1:nn
                offsets = x - gp.centres(:,i) * ones(1,nx);
                cov = gp.covariances(:,:,i);
                drive(i,:) = exp(-1/2 * sum(offsets .* (cov \ offsets), 1));
            end
        end
        
    end
    
end