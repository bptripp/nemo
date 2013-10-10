% Known issues:
% - the integration method for spiking simulations can cause delays of up
%   to one Bernoulli time
% - Gaussian correlations are corrected for spiking (based on 30Hz spiking 
%   for 100s) but not for rates
classdef PoissonSpikeGenerator < SpikeGenerator 
    
    properties (Access = private)
        intercepts = [];
        scaleFactors = [];
        cholCorrelations = [];
        cholSpikeCorrelations = [];

        spikeLastTime = [];
        spikeBernoulliTime = [];
        spikeA1 = [];
        spikeB1 = [];        
        spikeState1 = [];        
        spikeA2 = [];
        spikeB2 = [];
        spikeState2 = [];
    end
    
    properties (SetAccess = private)
        correlations = [];
        spikeCorrelations = [];
    end
    
    methods (Access = public)

        % correlations: correlations of underlying Gaussian distributions
        function psg = PoissonSpikeGenerator(intercepts, maxRates, correlations)
            psg = psg@SpikeGenerator(length(intercepts));

            assert(length(intercepts) == length(maxRates), 'Expected one intercept and one max rate per neuron')
            assert(size(intercepts, 1) == 1, 'Expected row vectors for maxRates, intercepts')
            assert(size(maxRates, 1) == 1, 'Expected row vectors for maxRates, intercepts')
            assert(size(correlations, 1) == size(correlations, 2), 'Correlation matrix should be square')
            assert(length(intercepts) == size(correlations, 1), 'Correlation matrix should have same dimension as # of intercepts, maxRates')
            
            psg.intercepts = intercepts;
            psg.scaleFactors = maxRates ./ (1 - intercepts);
            psg.correlations = correlations;
%             [V, D] = eig(correlations);
%             plot(diag(D)), pause
            psg.cholCorrelations = chol(correlations);
            
            psg.spikeCorrelations = correlations;
            psg.cholSpikeCorrelations = psg.cholCorrelations;
%             psg.spikeCorrelations = 2.1677*correlations + -1.0724*correlations.^2; % fit based on 30 Hz spiking 
%             psg.spikeCorrelations = atan(psg.spikeCorrelations);
%             psg.spikeCorrelations = max(-.99, min(.99, psg.spikeCorrelations));
            psg.spikeBernoulliTime = .001;

            Wn1 = 2*psg.spikeBernoulliTime*15; %15 Hz cutoff for correlated input
            [psg.spikeB1, psg.spikeA1] = butter(1, Wn1, 'low');
            
            Wn2 = 2*psg.spikeBernoulliTime*500;
            if Wn2 < 1
                [psg.spikeB2, psg.spikeA2] = butter(1, [Wn1 Wn2]);
            else 
                [psg.spikeB2, psg.spikeA2] = butter(1, Wn1, 'high');
            end
            
            psg.reset();
        end
        
        function setCorrelations(psg, correlations)
            psg.correlations = correlations;
            psg.cholCorrelations = chol(correlations);
            psg.spikeCorrelations = correlations;
            psg.cholSpikeCorrelations = psg.cholCorrelations;
            psg.reset();
        end

        % See SpikeGenerator.getRates(...). If startTime == endTime
        % then idealized rates are returned. If endTime > startTime
        % then the rates are based random samples from Poisson distributions
        % of counts for the interval endTime-startTime. 
        function rates = getRates(psg, drive, startTime, endTime)
            foo = ones(1, size(drive, 2));
            rates = max(0, (psg.scaleFactors'*foo) .* (drive - (psg.intercepts'*foo)));
            
            if (endTime > startTime)
                T = endTime - startTime;
                gaussianSamples = PoissonSpikeGenerator.randncov(size(rates, 2), psg.correlations, psg.cholCorrelations);
                cumProbs = normcdf(gaussianSamples, 0, 1);
                indices = find(rates > 0);
                rates(indices) = poissinv(cumProbs(indices), rates(indices)*T) / T;
            end
        end            
        
        % see SpikeGenerator.integrate(...)
        % Note that startTime is ignored except in the first call following a reset.
        function spikes = integrate(psg, drive, startTime, endTime)
            if ~isempty(psg.spikeLastTime)
                startTime = psg.spikeLastTime;
            end
            T = endTime - startTime;
            dt = psg.spikeBernoulliTime;
            n = floor(T/dt);
            psg.spikeLastTime = startTime + dt*n;
%             sprintf('n %i  start %f  end %f  last %f', n, startTime, endTime, psg.spikeLastTime)

            spikes = zeros(length(psg.intercepts), 1);
            if n > 0
                normSamples1 = PoissonSpikeGenerator.randncov(n, psg.spikeCorrelations, psg.cholSpikeCorrelations);
                [filtSamples1, psg.spikeState1] = filter(psg.spikeB1, psg.spikeA1, normSamples1, psg.spikeState1, 2);
                
                normSamples2 = randn(length(psg.scaleFactors), n);
                [filtSamples2, psg.spikeState2] = filter(psg.spikeB2, psg.spikeA2, normSamples2, psg.spikeState2, 2);

%                 filtSamples = filtSamples1 + filtSamples2;
                filtSamples = 3*filtSamples1 + filtSamples2;
                filtSamples = filtSamples / 1.165; %this corrects the rate based on above filtering and scaling
                
                rates = max(0, psg.scaleFactors' .* (drive - psg.intercepts'));
                p = rates * dt; % spike probabilities per time length dt
                cd = normcdf(filtSamples, 0, 1); % cumulative density
                spikes = cd > 1-p*ones(1,n); 
            end
            
            spikes = sum(spikes, 2);
        end
        
        function reset(psg)
            reset@SpikeGenerator(psg);
            psg.spikeLastTime = [];
            psg.spikeState1 = max(length(psg.spikeA1), length(psg.spikeB1)) - 1;
            psg.spikeState2 = max(length(psg.spikeA2), length(psg.spikeB2)) - 1;
        end
        
    end
    
    methods (Static) 

        % centres: values of encoded variable at which tuning curves peak
        % f: mean pairwise correlation as a function of two centres
        % sd: standard deviation of normal spread of correlations around f
        % correlations: resulting correlation matrix
        function correlations = genCorrelations(centres, f, sd)
            n = length(centres);
            lambda = 1;
            latentCentres = exprnd(lambda, n, 1);
%             latentCentres = rand(size(centres));
            correlations = zeros(n, n);
            for i = 1:n
                correlations(i,i) = 1;
                for j = (i+1):n
%                     cDist = norm(centres(:,i) - centres(:,j)); % distance between centres
                    lDist = abs(latentCentres(i)-latentCentres(j)); % distance between latent variables
%                     lDistCDF = 1-(2*lDist-lDist.^2); % cumulative density of lDist
                    lDistCDF = 1-exp(-lDist/lambda); % cumulative density of lDist                    
                    correlations(i,j) = norminv(1-lDistCDF , f(centres(:,i), centres(:,j)), sd); 
                    correlations(j,i) = correlations(i,j);
                end
            end

%             plotCorrelations(centres, latentCentres, correlations);
        end
        
        % Change eigenvalues to >0 (using Rebonato's method?) 
        % TODO: consider @article{Higham,
        % author = "Nicholas J. Higham",
        % journal = "IMA Journal of Numerical Analysis",
        % title = "Computing the nearest correlation matrix - a problem from finance",
        % year = 2002,
        % volume = 22,
        % pages = "329-343"
        % }
        function correlations = makePositiveDefinite(correlations)
            [V, D] = eig(correlations);
            
            minEig = 1e-10;
            D = diag(max(minEig, diag(D)));
            B1 = D.^(1/2) * V';
            
            %scale ... 
            for i = 1:size(B1, 2)
                B1(:,i) = B1(:,i) / norm(B1(:,i));
            end
            
            correlations = B1' * B1;
%             correlations = V * D * V';
        end

        % centres: tuning curve centres
        % correlations: matrix of spike count correlations
        % prob (optional): probability that a given pair will be plotted (to reduce
        %   plot density for large populations)
        % color (optional): plot color
        function plotCorrelationVsDistance(centres, correlations, varargin)            
            n = size(centres, 2);
            dim = size(centres, 1);
            points = (n^2 - n)/2;
            distances = zeros(dim, points);
            corrs = zeros(1, points);
            c = 0;
            for i = 1:n
                for j = (i+1):n
                    c = c + 1;
                    distances(:,c) = abs(centres(:,i) - centres(:,j));
                    corrs(c) = correlations(i,j);
                end
            end
            
            prob = 1;
            newRandStream = RandStream('mt19937ar', 'Seed', 2000);
            defaultRandStream = RandStream.setDefaultStream(newRandStream);
            if ~isempty(varargin)
                prob = varargin{1};
            end
            ind = find(rand(size(distances,2), 1) <= prob);
            RandStream.setDefaultStream(defaultRandStream);
            
            colour = 'b';
            if length(varargin) > 1
                colour = varargin{2};
            end

            if dim == 1
                scatter(distances(ind), corrs(ind), [], colour);
                xlabel('Pairwise Difference Between Tuning Centres')
                ylabel('Pairwise Correlation')
            elseif dim == 2
                scatter3(distances(1,ind), distances(2,ind), corrs(ind), [], colour);
                xlabel('Pairwise Difference Between Tuning Centres')
                zlabel('Pairwise Correlation')
            elseif dim == 3
                subplot(1, 3, 1)
                scatter3(distances(1,ind), distances(2,ind), corrs(ind), [], colour);
                zlabel('Pairwise Correlation')
                subplot(1, 3, 2)
                scatter3(distances(1,ind), distances(3,ind), corrs(ind), [], colour);
                subplot(1, 3, 3)
                scatter3(distances(2,ind), distances(3,ind), corrs(ind), [], colour);
            else 
                error('Can''t plot correlations for more than 3 represented dimensions')
            end
        end
        
        % Returns samples from a multivariate normal distribution with given
        % covariance. 
        % 
        % n: number of (vector) samples
        % COV: covariance matrix
        % cholCOV (optional): Cholesky decomposition of COV, which can be
        %   provided to save time if this method is called repeatedly with
        %   the same covariance matrix
        % result: n random samples 
        % R (optional output): correlation matrix of samples
        function [result, varargout] = randncov(n, COV, varargin)
            %TODO: deal with correlations of 1 (maybe cap correlation and copy afterward)
            [a, b] = size(COV);
            
            if ~isempty(varargin) 
                cholCOV = varargin{1};
            else 
                cholCOV = chol(COV);
            end
            
            result = cholCOV' * randn(a, n);
            
            if nargout > 1 
                varargout{1} = corrcoef(result');
            end        
        end

    end
end

% Plots correlations generated using a latent variable to show how they 
% vary with distance between tuning curve centres and the latent variable 
% -- see genCorrelations(...) 
function plotCorrelations(centres, latentCentres, correlations) 
    n = length(latentCentres);
    correctedCorr = CorrCountModel.makePositiveDefinite(correlations);
    max(max(abs(correctedCorr - correlations)))
    num = (n^2 - n) / 2;
    ld = zeros(1, num);
    cd = zeros(1, num);
    cc = zeros(1, num);
    ccc = zeros(1, num);

    c = 0;
    for i = 1:n
        for j = (i+1):n
            c = c + 1;
            ld(c) = norm(latentCentres(i) - latentCentres(j));
            cd(c) = norm(centres(:,i) - centres(:,j));
            cc(c) = correlations(i,j);
            ccc(c) = correctedCorr(i,j);
        end
    end

    figure, hold on
    scatter3(ld, cd, cc, 'b');
    scatter3(ld, cd, ccc, 'r');
end
        
