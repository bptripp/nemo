% Plots tuning curves of a population along a single axis. 
% 
% population: Population to plot
% x (optional): Vector value at which to plot the tuning curves. NaN
%   elements are taken as dimensions along which to plot the curves, with
%   real-valued elements indicating constants in other dimensions. 
%   At least one element must be NaN (no more than 3). Default is all NaNs. 
% T (optional): integration time (for intrinsically noisy populations;
%   defaults to zero)
% indices (optional): neuron indices (defaults to all for 1D and the first
%   nine for >1D)
% plotColor (optional): plot color
function plotPopulation(population, varargin)
    x = NemoUtils.getOptionalArg(varargin, 1, NaN, 'x', [population.dimension 1]);
    plotColor = NemoUtils.getOptionalArg(varargin, 4, '', 'plotColor', []);
%     x = NaN * ones(population.dimension, 1);
%     if ~isempty(varargin) && ~isempty(varargin{1})
%         x = varargin{1};
%         assert(size(x, 2) == 1, 'Expected column vector for x');
%         assert(size(x, 1) == length(population.radii), 'x must have same dimension as population');
%     end
    nPlotDims = sum(isnan(x));

    if nPlotDims == 1
        indices = 1:population.spikeGenerator.n;
    else 
        indices = 1:min(population.spikeGenerator.n, 9);
    end
    
    T = NemoUtils.getOptionalArg(varargin, 2, 0, 'T', 1);
    indices = NemoUtils.getOptionalArg(varargin, 3, indices, 'indices', []);
    indices = indices(indices <= population.spikeGenerator.n);
%     T = 0;
%     if ~isempty(varargin)
%         if length(varargin) > 1 && ~isempty(varargin{2})
%             T = varargin{2};
%         end
%         if length(varargin) > 2 && ~isempty(varargin{3})
%             indices = varargin{3};
%             indices = indices(indices <= population.spikeGenerator.n);
%         end
%     end

    if nPlotDims == 1
        plotPopulation1D(population, x, T, indices, plotColor);
    elseif nPlotDims == 2
        plotPopulation2D(population, x, T, indices, plotColor);  
    elseif nPlotDims == 3
        plotPopulation3D(population, x, T, indices, plotColor); 
    else 
        error('Can''t plot less than one or more than two dimensions')
    end
end

function plotPopulation1D(population, x, T, indices, plotColor) 
    dim = find(isnan(x));
    r = population.radii(dim);    
    points = population.offsets(dim) + (-r:(r/50):r);
    x = x * ones(1, length(points));
    x(dim,:) = points;
        
    rates = getRates(population, x, 0, T);
    
    figure
    plot(points, rates(indices,:)', plotColor)
    title(sprintf('Population %s dimension %i', population.name, dim));
    xlabel('represented variable')
    ylabel('firing rate')
end

function plotPopulation2D(population, x, T, indices, plotColor)
    dims = find(isnan(x));
    r1 = population.radii(dims(1));
    points1 = -r1:(r1/20):r1;
    n1 = length(points1);
    r2 = population.radii(dims(2));
    points2 = -r2:(r2/25):r2;
    n2 = length(points2);

    points = [reshape(points1' * ones(1, n2), 1, n1*n2); reshape(ones(n1, 1) * points2, 1, n1*n2)]; 
    x = x * ones(1, size(points, 2));
    x(dims,:) = points;
    rates = getRates(population, x, 0, T);
    figure
    side = ceil(length(indices)^.5);
    for i = 1:length(indices)
        subplot(side, side, i);
        mesh(points1, points2, reshape(rates(indices(i), :), n1, n2)');
    end
end

function plotPopulation3D(population, x, T, indices, plotColor)
    dims = find(isnan(x));
    r1 = population.radii(dims(1));
    points1 = -r1:(r1/5):r1;
    r2 = population.radii(dims(2));
    points2 = -r2:(r2/5):r2;
    r3 = population.radii(dims(3));
    points3 = -r3:(r3/5):r3;

    [X, Y, Z] = meshgrid(points1, points2, points3);
    points = [reshape(X, 1, numel(X)); reshape(Y, 1, numel(Y)); reshape(Z, 1, numel(Z))];

    x = x * ones(1, size(points, 2));
    x(dims,:) = points;
    rates = getRates(population, x, 0, T);
    sizes = 1 + 143 * rates / max(max(rates));
    figure
    side = ceil(length(indices)^.5);
    for i = 1:length(indices)
        subplot(side, side, i);
        scatter3(points(1,:), points(2,:), points(3,:), sizes(indices(i), :), sizes(indices(i), :))
    end
end

