% A population output that is decoded from spiking activity. An estimate of 
% some function of the variable represented by the population. 
classdef DecodedOrigin < Origin
    
    properties (SetAccess = private)
        population = [];
        f = [];
    end
    
    properties (SetAccess = protected)
        decoders = [];
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
        
        % decoders: Linear decoding vectors for approximation of some function of the variable 
        %   represented by the associated population activity
        function setDecoders(do, decoders)
            assert(size(decoders, 1) == do.dim, 'Each column of decoders should have length %i (not %i)', do.dim, size(decoders,1)); 
            do.decoders = decoders;
        end
        
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
%             relNoise = NemoUtils.getOptionalArg(varargin, 4, .2, 'relNoise', 1);
            relNoise = NemoUtils.getOptionalArg(varargin, 4, .05, 'relNoise', 1);
            assert(maxCorrChange < 1, 'maxCorrChange should be between 0 and 1')
            
            if ~isempty(varargin) && ~isempty(varargin{1})
                points = varargin{1};
                assert(size(points, 1) == length(do.population.radii), 'Expected points of same dimension as population');
            else 
                pr = do.population.radii;
                n = 300 + 1200*(length(pr)>1);
                points = Population.genRandomPoints(n, pr, 0, do.population.ellipsoidRegion, do.population.offsets);
%                 points = Population.genRandomPoints(n, pr, 2, do.population.ellipsoidRegion, do.population.offsets);
            end

%             d = length(population.radii);
%             [n, points] = DecodedOrigin.getPoints(d, varargin{:});          
%             figure, scatter(points(1,:), points(2,:)), title(sprintf('%i points', n)), pause
            
            ideal = do.f(points);  
            rates = do.getNoisyRates(points, T, 0);            
            gamma = rates * rates';
            gamma = gamma + (relNoise*max(rates(:)))^2*size(rates,2) * eye(size(gamma,1));
            
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
            
            if nargout > 0
                varargout{1} = norms;
            end
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
        function setActivity(do, time, activity)
            assert(~isempty(do.decoders), 'Either setDecoders(...) or findDecoders(...) must be called before this DecodedOrigin is used');
            assert(size(activity, 2) == 1, 'Activity should be a column vector');
            assert(size(activity, 1) == size(do.decoders, 2), 'Activity should be a vector of length %i, not %i', ...
                size(do.decoders, 2), size(activity, 1));

            output = do.decoders * activity;
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
