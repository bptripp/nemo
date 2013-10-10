% Plots the static decoding error in at most two dimensions of a
% DecodedOrigin. 
% 
% population: A Population
% originName: The name of the DecodedOrigin on the above population that is
%   to be plotted
% dims (optional): a list of one or two axes along which to plot the error, 
%   if the decoded function is more than two-dimensional 
% outDim (optional): dimension of output to plot (one at a time; default 1) 
% indices (optional): indices of active neurons and/or clusters (default 1:n)
function mse = plotDecodedOrigin(population, originName, varargin)
    if length(population.radii) == 1
        dims = 1;
    else 
        dims = [1; 2];
    end
    if ~isempty(varargin) && ~isempty(varargin{1})
        dims = varargin{1};
        assert(length(dims) == 1 || length(dims) == 2, 'Dimensions to plot should have length 1 or 2, not %i', length(dims));
    end
    indices = NemoUtils.getOptionalArg(varargin, 3, 1:population.spikeGenerator.n, 'indices', []);
    
    origin = [];
    for i = 1:length(population.origins)
        if strcmp(originName, population.origins{i}.name) 
            origin = population.origins{i};
        end
    end
    if isempty(origin)
        for i = 1:length(population.origins) 
            sprintf('Origin: %s', population.origins{i}.name)
        end
        error('Origin %s not found', originName);
    end

    outDim = NemoUtils.getOptionalArg(varargin, 2, 1, 'outDim', 1);
    
    if length(dims) == 1
        r = population.radii(dims);
        points = population.offsets(dims) + (-r:(r/50):r);
        ideal = origin.f(points);
        ideal = ideal(outDim,:);
        x = zeros(length(population.radii), length(points));
        x(dims,:) = points;
        approx = getOriginOutput(origin, x, indices);
        approx = approx(outDim,:);
%         rates = getRates(population, x, 0, 0);
%         approx = origin.decoders * rates; 
        mse = mean((approx-ideal).^2)^.5;
        figure, hold on
        plot(points, ideal, 'k')
        plot(points, approx, 'k--')
        xlabel('represented variable')   
        ylabel('decoded function')
        ha2 = axes('XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XTick', [], 'YColor', 'r', 'Position', get(gca, 'Position'));   
        line(points, approx - ideal, 'Color', 'r', 'Parent', ha2)
        ylabel('decoding error')        
    elseif length(dims) == 2
        r = population.radii(dims);
        x1 = -r(1):(r(1)/20):r(1);
        x2 = -r(2):(r(2)/30):r(2);
        [X1, X2] = meshgrid(x1, x2);
        points = [reshape(X1, 1, numel(X1)); reshape(X2, 1, numel(X2))];
        x = zeros(length(population.radii), size(points,2));
        x(dims,:) = points;
        ideal = origin.f(x);
        ideal = ideal(outDim,:);
        approx = getOriginOutput(origin, x, indices);
        approx = approx(outDim,:);
        mse = mean((approx-ideal).^2)^.5;
%         rates = getRates(population, x, 0, 0);
%         approx = origin.decoders * rates;
        figure, mesh(X1, X2, reshape(ideal, length(x2), length(x1))), zlabel('ideal')
        figure, mesh(X1, X2, reshape(approx, length(x2), length(x1))), zlabel('decoded')
        figure, mesh(X1, X2, reshape(approx-ideal, length(x2), length(x1))), zlabel('error'), title(sprintf('RMS %f', mse))
    else 
        error('Can''t plot decoding over >2 dimensions')
    end
    
end

function output = getOriginOutput(o, x, indices)
    rates = getRates(o.population, x, 0, 0, indices);
    nx = size(rates, 2);
    output = zeros(o.dim, nx);
    for i = 1:nx
        o.setActivity(0, rates(:,i), indices);
        output(:,i) = getOutput(o);
    end
end
