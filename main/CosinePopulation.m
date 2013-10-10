% A population of neurons with cosine tuning.  
classdef CosinePopulation < Population
    properties (Access = public)
        encoders = [];
    end
    
    properties (SetAccess = private)
%         Z = [];
%         T = [];
%         clusterRadialDrive = [];
%         clusterRFs = [];
%         clusterEncoders = [];
%         clusterCosts = [];
%         clusterIndices = [];
% %         clusterErrors = [];
%         originErrors = struct();
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
        function drive = getDrive(cp, x, varargin)
            if ~isempty(cp.offsets)
                x = x - repmat(cp.offsets, 1, size(x, 2));
            end
            
            if ~isempty(varargin) % indices are specified   
                drive = cp.encoders(:,varargin{1})' * x;
            else 
                drive = cp.encoders' * x;
            end
        end
        
        % TODO: abstract method for Population
        function initClusters(cp)
            n = cp.spikeGenerator.n;
            nPoints = 3000;
            points = Population.genRandomPoints(nPoints, cp.radii, 0);
            rates = getRates(cp, points, 0, 0);
            rho = corr(rates');
            rho = rho + diag(NaN*ones(1, n));
            rho = [rho NaN*ones(n,n-1); NaN*ones(n-1,n+n-1)];
            
            available = 1:n; %available to merge
            cp.Z = zeros(n-1, 2); %TODO: need this information to nagivate tree, but maybe a better way to store it
            rates = [rates; zeros(n-1, nPoints)];
            for i = 1:n-1
                i
                ind = find(rho == max(rho(:)), 1, 'first');
                row = ceil(ind/size(rho, 1));
                col = ind - (row-1)*size(rho, 1); 
                cp.Z(i,:) = [row col];
                addMerge(cp, [row col]);
                
                ratesi = getRates(cp, points, 0, 0, n+i);
                rates(n+i,:) = ratesi;
                
                available = [setdiff(available, [row col])];
                rhoi = corr(ratesi', rates(1:n+i-1,:)');
%                 rho(n+i,1:n+i-1) = rhoi;
                rho(n+i,available) = rhoi(available);
                available = [available n+i];
                
                % merged clusters are not considered for future merges ...
                rho(row,:) = NaN;
                rho(:,row) = NaN;
                rho(:,col) = NaN;
                rho(col,:) = NaN;
            end
        end
        
        % Creates a merged neuron with properties that are intermediate to
        % those of the given indices. This is to support reduced
        % simulations. 
        % 
        % indices: indices of neurons from which to create a new neuron
        %   that has properties that are roughly the average of those in
        %   the merged group. Indices <=n are neurons, and those >n are
        %   past merges (in which case the property averages should be
        %   weighted according to the size of previously merged groups). 
        function addMerge(cp, indices)
            % create group in spike generator
            addMerge(cp.spikeGenerator, indices);

            % create new encoder from weighted average
            counts = cp.spikeGenerator.mergedCount(indices);
            
            % stretch existing encoders to hypersphere, average, normalize,
            % un-stretch to 1/radii
            encodersToMerge = cp.encoders(:,indices);
            encodersToMerge = encodersToMerge .* repmat(cp.radii, 1, length(indices));
            newEncoder = encodersToMerge * counts / sum(counts);            
            newEncoder = newEncoder ./ cp.radii / norm(newEncoder);
            cp.encoders = [cp.encoders newEncoder];

            % sum decoders in all origins
            for i = 1:length(cp.origins)
                if isa(cp.origins{i}, 'DecodedOrigin')
                    addMerge(cp.origins{i}, indices);                
                end
            end
            
            % TODO: calculate error increment in origins due to merge
        end
        
        function removeMerge(cp, indices)
            assert(min(indices) > cp.spikeGenerator.n, 'Can only use this method to remove merge neurons');

            removeMerge(cp.spikeGenerator, indices);

            toKeep = setdiff(1:cp.spikeGenerator.n, indices);
            cp.encoders = cp.encoders(:,toKeep);

            for i = 1:length(cp.origins)
                if isa(cp.origins{i}, 'DecodedOrigin')
                    removeMerge(cp.origins{i}, indices);                
                end
            end
        end
    
%         function initClusters(cp)
%             n = cp.spikeGenerator.n;
%             cp.Z = linkage(cp.encoders');
%             
%             % create clustered response function for each merge
%             cp.clusterRadialDrive = [-2 -1.5 -1:.1:1 1.5 2];
%             RF = zeros(2*cp.spikeGenerator.n-1, length(cp.clusterRadialDrive));
%             
%             for i = 1:length(cp.clusterRadialDrive)
%                 RF(1:n,i) = run(cp.spikeGenerator, repmat(cp.clusterRadialDrive, n, 1), 0, 0, 0);
%             end
%             
%             cp.clusterIndices = cell(n-1, 1);
%             for i = 1:(size(cp.Z, 1)-1) 
%                 ind1 = cp.Z(i,1);
%                 ind2 = cp.Z(i,2);
%                 RF(n+i,:) = RF(ind1,:) + RF(ind2,:); 
%                 cp.clusterIndices{i} = [cp.clusterIndices{ind1} cp.clusterIndices{ind2}];
%                 cp.clusterEncoders = mean(cp.encoders(:,cp.clusterIndices{i}), 2);
%             end
%             cp.clusterRFs = RF(n+1:end,:);            
%             
%             cp.T = cluster(Z, 'MaxClust', 1:n);
%         end

        %TODO: test
        function initOriginErrors(cp)
            points = Population.genRandomPoints(1500, cp.radii, 0);
            rates = getRates(cp, points, 0, 0, 1:cp.spikeGenerator.n+size(cp.Z,1));
            for i = 1:length(cp.origins)
                if isa(cp.origins{i}, 'DecodedOrigin')
                    findErrorEstimates(cp.origins{i}, rates, cp.Z);
%                     initOriginError(cp, cp.origins{i}.name);
                end
            end
        end
        
%         % originName: 
%         function initOriginError(cp, originName)
%             %TODO: do I have to simulate all cluster combinations to model noise?
%             origin = getOrigin(cp, originName);
%             points = Population.genRandomPoints(1500, cp.radii, 0);
%             
%             setfield(cp.originErrors, originName, errors);
%         end
%         
%         % tolerances: struct with origin names and vectors of error 
%         %   tolerances per dimension of corresponding origins
%         function indices = chooseClusters(cp, tolerances)
%             
%         end

    end
end
