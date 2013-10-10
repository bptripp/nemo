% A population output that is decoded from spiking activity. An estimate of 
% some function of the variable represented by the population. 

% TODO: variable dimension; error in synaptic current as a function of dimension

classdef DecodedOrigin < Origin
    
    properties (SetAccess = private)
        population = [];
        f = [];
        PCDecoders = [];        
    end
    
    properties (SetAccess = protected)
        decoders = [];
        tolerance = [];
        mergeErrors = []; %error variance associated with each merge, relative to no merges (n-1 by dim)
        splitValues = []; %reduction in error with split, divided by # splits needed to achieve it (error may not drop 
                          %   monotonically, so we keep splitting in order
                          %   of clusters (order of signal correlations)
                          %   until error is reduced in each dimension, and
                          %   take the best value along this path
    end
    
    properties (Access = private)
        activityU = [];
        activityS = [];
    end
    
    methods (Access = public)

        % Either setDecoders(...) or findDecoders(...) must be called
        % before the object is used. 
        % 
        % name: A name that distinguishes this origin from others of the same
        %       Node, e.g. different decoded outputs. 
        % f: A function of the population's represented domain that is to be decoded
        %       by the origin
        % population: Population to which the origin belongs
        function do = DecodedOrigin(name, f, population)
            dim = length(population.radii);
            assert(strcmp(class(f), 'function_handle'), 'Expected a function handle');
            testOutput = f(zeros(dim, 100));
            assert(size(testOutput, 2) == 100, 'Function output should be a matrix with one column vector for each input (multiple inputs may be given at once in separate columns)')
            outputDim = size(testOutput, 1);
                        
            do = do@Origin(name, outputDim);
            do.f = f;
            do.population = population;
        end
        
        % Finds the neurons that are important for determining selected
        % dimensions of the output. 
        %
        % dimensions: list of dimensions of interest
        % tolerance: fraction of decoder ??? that can be ignored
        % indices: indices of neurons that are important contributors to
        %   the given dimensions
%         function indices = getRelevantNeurons(dimensions, tolerance)
%             %TODO: implement; choose efficient thresholding method
%         end
        
        % decoders: Linear decoding vectors for approximation of some function of the variable 
        %   represented by the associated population activity
        function setDecoders(do, decoders)
            assert(size(decoders, 1) == do.dim, 'Each column of decoders should have length %i (not %i)', do.dim, size(decoders,1)); 
            do.decoders = decoders;
            updatePCDecoders(do);            
        end
        
        function setTolerance(do, tolerance) 
            assert(length(tolerance) == size(do.decoders, 1), 'Expected one tolerance per output dimension');
            do.tolerance = tolerance;
        end
        
        % To be called by the population after clusters have been set up.
        % Calculates error variance associated with each merge, if no other
        % neurons are merged. Error is relative to un-merged output rather
        % than the ideal do.f(...). 
        % 
        % Sets fields mergeErrors and splitValues. 
        % 
        % TODO: looks right but how can I unit-test this? simple known data?
        %
        % rates: rates of neurons and clusters at sample points
        % Z: cluster hierarchy
        function findErrorEstimates(do, rates, Z)
            n = do.population.spikeGenerator.n;
            dim = size(do.decoders, 1);
            assert(size(Z, 1) == n-1, 'Expected Z to have n-1 rows')
            ideal = do.decoders(:,1:n) * rates(1:n,:);
            
            do.mergeErrors = zeros(dim, n-1);
            for i = 1:n-1
                activeIndices = [n+i setdiff(1:n, DecodedOrigin.getChildren(n+i, Z))];
                approx = do.decoders(:,activeIndices) * rates(activeIndices, :);
                do.mergeErrors(:,i) = mean((approx - ideal).^2, 2); %note: mean is not typically zero
            end
            
            % update split value function ... 
            % TODO: normalize variance or SD by cost?
            do.splitValues = zeros(dim, n-1); % splitImprovements ./ splitCosts
            splitImprovements = zeros(size(do.splitValues));
            splitCosts = zeros(size(do.splitValues)); % number of splits needed to obtain above improvements
            for i = 1:n-1
                children = Z(i,:);
                improvements = zeros(dim, 3); % col 1 = child 1; col 2 = child 2; col 3 = immediate
                costs = ones(dim, 3);
                improvements(:,3) = do.mergeErrors(:,i);
                if children(1) > n
                    improvements(:,1) = splitImprovements(:,children(1)-n);
                    improvements(:,3) = improvements(:,3) - do.mergeErrors(:,children(1)-n);
                    costs(:,1) = costs(:,1) + splitCosts(:,children(1)-n);
                end
                if children(2) > n
                    improvements(:,2) = splitImprovements(:,children(2)-n);
                    improvements(:,3) = improvements(:,3) - do.mergeErrors(:,children(2)-n);
                    costs(:,2) = costs(:,2) + splitCosts(:,children(2)-n);
                end
                
                values = improvements ./ costs;
                
                [Y, maxInd] = max(values, [], 2);
                for j = 1:dim
                    splitImprovements(j,i) = improvements(j,maxInd(j));
                    splitCosts(j,i) = costs(j,maxInd(j));
                end
                do.splitValues(:,i) = splitImprovements(:,i) ./ splitCosts(:,i);
            end
        end
        
        % indices: indices of neurons and clusters to be simulated together
        %   (should exactly span all neurons)
        % error: error in each output dimension associated with given set
        %   of clusters
        function error = getErrorEstimate(do, indices)
            n = do.population.spikeGenerator.n;
            clusterIndices = indices(indices > n) - n;
            error = sum(do.mergeErrors(:,clusterIndices), 2);
        end
        
%         %TODO: docs
%         function setClusterDecoders(do, Z)
%             assert(size(Z, 1) == size(do.decoders, 2)-1, sprintf('Expected %i merges', size(do.decoders, 2)-1));
%             %TODO: this is the wrong approach ... have to store RFs
%             %multiplies by decoders here
%             do.clusterDecoders = zeros(size(do.decoders, 1), size(do.decoders, 2)-1); % #dims X #merges
%             for i = 1:size(Z, 1)
%                 ind1 = Z(i,1);
%                 ind2 = Z(i,2);
%             end
%         end
        
        % points (optional): points at which firing rates are sampled
        %   in order to find decoders (defaults to a uniform distribution
        %   in the represented domain, the number of which increases with
        %   dimension somewhat more slowly than the volume). Reasons for
        %   specifying these points include wanting more points than the
        %   default, and wanting a different density in order to emphasize
        %   accuracy in certain regions. 
        % maxCorrChange (optional): if >0, rates will be obtained
        %   repeatedly for each point until the largest change in pairwise 
        %   rate correlation from one iteration to the next is 
        %   less than maxCorrChange
        % T (optional): run time (default 0)
        % relNoise (optional): noise added to rates to make optimal
        %   decoders less sensitive (proportion of highest rate; default .05)
        % 
        % dCorr (optional output): history of largest correlation changes
        function varargout = findDecoders(do, varargin)
            maxCorrChange = NemoUtils.getOptionalArg(varargin, 2, 0, 'maxCorrChange', 1);
            T = NemoUtils.getOptionalArg(varargin, 3, 0, 'T', 1);
            relNoise = NemoUtils.getOptionalArg(varargin, 4, .05, 'relNoise', 1);
            assert(maxCorrChange < 1, 'maxCorrChange should be between 0 and 1')
            
            if ~isempty(varargin) && ~isempty(varargin{1})
                points = varargin{1};
                assert(size(points, 1) == length(do.population.radii), 'Expected points of same dimension as population');
            else 
                pr = do.population.radii;
                n = 300 + 1200*(length(pr)>1);
                points = Population.genRandomPoints(n, pr, 0, do.population.ellipsoidRegion, do.population.offsets);
            end

%             d = length(population.radii);
%             [n, points] = DecodedOrigin.getPoints(d, varargin{:});          
%             figure, scatter(points(1,:), points(2,:)), title(sprintf('%i points', n)), pause
            
            ideal = do.f(points);  
            rates = do.getNoisyRates(points, T, relNoise);
            gamma = rates * rates';
            dCorr = [];
            V = rates * ideal';
            i = 1;
%             for i = 2:10
            while maxCorrChange > 0 && (isempty(dCorr) || dCorr(end) > maxCorrChange)
                i = i + 1;
                rates = do.getNoisyRates(points, T, relNoise);
                newGamma = gamma + rates * rates';
                dCorr = [dCorr maxCorrDiff(newGamma, gamma)];
                gamma = newGamma;
                V = V + rates * ideal';
                sprintf('Iteration %i: Largest correlation change %f', 1+length(dCorr), dCorr(end))
            end
            invgamma = pinv(gamma);
            do.decoders = (invgamma*V)';
            
            updatePCDecoders(do);
            
            if nargout > 0
                varargout{1} = norms;
            end
        end

        % U: first part of singular value decomposition of activities (neuron by sample)
        % S: second part of SVD
        function initPC(do, U, S)
            assert(~xor(isempty(U), isempty(S)));
            do.activityU = U;
            do.activityS = S;
            updatePCDecoders(do);
        end
        
        % To be called when neurons are merged in support of reduced
        % simulations. 
        % 
        % indices: indices of neurons that have been merged in the
        %   population. Indices <=n are neurons, and those >n are
        %   past merges (in which case the property averages should be
        %   weighted according to the size of previously merged groups). 
        function addMerge(do, indices)
            do.decoders = [do.decoders sum(do.decoders(:,indices), 2)];
        end
        
        % Removes a merged neuron. 
        % 
        % index: indices of merges to remove (must be >n) 
        function removeMerge(do, indices)
            assert(min(indices) > do.population.spikeGenerator.n, 'Cannot use this method to remove actual neurons -- only merge neurons that model groups of actual neurons');
            toKeep = setdiff(1:size(do.decoders, 2), indices);
            do.decoders = do.decoders(:,toKeep);
        end

        function rates = getNoisyRates(do, points, T, relNoise)
            rates = getRates(do.population, points, 0, T);
            if (relNoise > 0)
                noise = max(max(abs(rates))) * relNoise * randn(size(rates));
                rates = rates + noise;
            end
        end
        
        % This should normally be called by the associated Node before 
        % completion of run(...), in order to make Node output available 
        % to other Nodes. When running in DIRECT_MODE, there is no activity
        % and setX(...) should be called instead. 
        % 
        % time: Simulation time (at end of run step)
        % activity: Firing rates of neurons in the population, or spikes with integral 1 (o.dim by 1)
        % indices (optional): indices of neurons and/or clusters for which
        %   activity is provided (default 1:n)
        function setActivity(do, time, activity, varargin)
            indices = NemoUtils.getOptionalArg(varargin, 1, 1:do.population.spikeGenerator.n, 'indices', []);
            
            assert(~isempty(do.decoders), 'Either setDecoders(...) or findDecoders(...) must be called before this DecodedOrigin is used');
            assert(size(activity, 2) == 1, 'Activity should be a column vector');
            assert(size(activity, 1) == length(indices), 'Activity should be a vector of length %i, not %i', ...
                length(indices), size(activity, 1));
            assert(max(indices) <= size(do.decoders, 2), 'There are only %i decoders', size(do.decoders, 2));

            output = do.decoders(:,indices) * activity;
%             size(output)
            setOutput(do, time, output);
        end

        % To be called by the associated Population before the end of a run(...)
        % when operating in DIRECT_MODE. 
        % 
        % time: Simulation time (at end of run step)
        % x: State variable represented by the associated Population
        function setX(do, time, x)
            output = do.f(x);
            setOutput(do, time, output);
        end
        
        % To be called by the associated Population before the end of a run(...)
        % when operating in PC_MODE. 
        % 
        % time: Simulation time (at end of run step)
        % x: Principal components of activities of the associated Population
        function setPC(do, time, PC)
            assert(~isempty(do.PCDecoders), 'Both initPC and either setDecoders or findDecoders must be called before this DecodedOrigin is used in PC mode');
            assert(size(PC, 2) == 1, 'PCs should be a column vector');
            assert(size(PC, 1) <= size(do.PCDecoders, 2), 'PCs should be a vector of length <= %i', size(do.PCDecoders, 2));
            
            output = do.PCDecoders * activity;
            setOutput(do, time, output);
        end
        
    end
    
    methods (Access = private)
        function updatePCDecoders(do) 
            if ~isempty(do.decoders) && ~isempty(do.acitivityU) && ~isempty(do.acitivityS)
                do.PCDecoders = do.decoders * do.acitivityU * do.acitivityS;
            end
        end
    end
    
    methods (Static) 

%         function [n, points] = getPoints(d, varargin)
%             if ~isempty(varargin) && ~isempty(varargin{1})
%                 points = varargin{1};
%                 n = size(points, 2);  
%                 assert(size(points, 1) == d, 'Expected points of same dimension as population');
%             else 
%                 n = 200^(2/3) * 200.^(d/3)
%                 points = Population.genRandomPoints(n, population.radii, false, population.ellipsoidRegion);
%             end
%         end
        
        % activity: firing rates of neurons in a population at selected
        %   points
        % ideal: values of the function to be decoded at the same points
        % relNoise: amount of noise to be added to activity, as a
        %   proportion of max activity (0 to 1)
        % decoders: optimal linear decoders (activity * decoders optimally
        %   approaches ideal)
%         function decoders = optimalDecoders(activity, ideal, relNoise)
%             if (relNoise > 0)
%                 noise = max(max(abs(activity))) * relNoise * randn(size(activity));
%                 activity = activity + noise;
%             end
% 
%             gamma = activity * activity'; 
%             invgamma = pinv(gamma);
% 
%             V = activity * ideal';
%             decoders = (invgamma*V)';
%         end
        
        %TODO: this is probablty slower as a static method
        function children = getChildren(cluster, Z)
            n = size(Z, 1) + 1;
            children = Z(cluster - n,:);
            if children(1) > n
                children = [children DecodedOrigin.getChildren(children(1), Z)];
            end
            if children(2) > n
                children = [children DecodedOrigin.getChildren(children(2), Z)];
            end    
        end        
        
        
    end
end

function result = maxCorrDiff(cov1, cov2)
    %TODO: these variable names are confusing (not really covariance and
    %SD but products with n plus mean^2)
    sd1 = diag(cov1).^(1/2);
    sd2 = diag(cov2).^(1/2);
    corr1 = cov1 ./ (sd1 * sd1');
    corr2 = cov2 ./ (sd2 * sd2');
    
    result = max(max(abs(corr1-corr2)));
end


% function children = getChildren(cluster, Z)
%     n = size(Z, 1) + 1;
%     children = Z(cluster - n,:);
%     if children(1) > n
%         children = [children getChildren(children(1), Z)];
%     end
%     if children(2) > n
%         children = [children getChildren(children(2), Z)];
%     end    
% end