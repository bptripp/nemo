% probe: A Probe that has collected state data (e.g. decoded output,
%   spikes)
% tau (optional): Time constant (s) for filtering output. The default is 
%   no filter (this can also be specified with an empty array []). Decoded 
%   spike output is typically filtered with a time constant of 0.005 or 
%   0.01, to approximate the current induced by the signal at fast synapses. 
% fh (optional): Figure handle, if it is desired to plot in a specified
%   figure. By default a new figure is created. 
% clf (optional): 1=clear current figure before plotting; 0=don't (default)
% color (optional): color spec, passed to Matlab line function (default 'b';
%   ignored for spikes)
% h: line handle (to allow caller to set properties)
function h = plotProbe(probe, varargin)
    if ~isempty(varargin) && ~isempty(varargin{1})
        tau = varargin{1};
        [time, history] = getHistory(probe, tau);
    else 
        [time, history] = getHistory(probe);
    end
    
    if length(varargin) < 2 || isempty(varargin{2})
        fh = figure;
    else 
        fh = figure(varargin{2});
    end
    
    if length(varargin) >= 3 && ~isempty(varargin{3}) && varargin{3}
        clf(fh);
    end
    
    spec = 'b';
    if length(varargin) >= 4 
        spec = varargin{4};
    end
        
    hold on
    if islogical(history) % plot spikes
        for i = 1:size(history, 1)
            spikeIndices = find(history(i,:));
            st = time(spikeIndices);
%             nn = ones(size(spikeTimes))*i;
% %             plot([spikeTimes; spikeTimes], [nn-1; nn+1], 'k')
%             plot(spikeTimes, nn, 'k.')
            h = line(st, i*ones(size(st)));
            set(h, 'LineStyle', 'none')
            set(h, 'Marker', 's')
            set(h, 'MarkerEdgeColor', 'k')
            set(h, 'MarkerFaceColor', 'k')
            set(h, 'MarkerSize', 1)
        end
        xlabel('spike times (s)')
        ylabel('neuron #')
        set(gca, 'YLim', [-.5 size(history, 1)+.5])
    else % plot vector
        h = line(time, history);
        set(h, 'Color', spec)
        xlabel('Time (s)', 'FontSize', 18)
        ylabel(probe.stateName, 'FontSize', 18)
        set(gca, 'FontSize', 18)
    end
end