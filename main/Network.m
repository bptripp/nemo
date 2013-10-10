% A collection of Nodes and connections between them, which run with a 
% common time step and exchange information at the end of each step. 
% Nodes in a Network run independently for steps of fixed length, and 
% exchange information at the beginning of each step.  
classdef Network < Node
    
    properties (Access = public)
        nodes = cell(0, 0); 
        connections = [];
        probes = [];
        step = []; %step size (s)
    end
    
    methods (Access = public)
        
        % step: constant step size (s) 
        function n = Network(step)
            n.step = step;
        end

        % start: beginning of simulation time (s)
        % stop: end of simulation time (s)
        function run(n, start, stop)
            %TODO: in reduced mode, may have to run specific state
            %variables (neurons with substantial decoders related to
            %specific outputs) multiple times with different fidelities
            %TODO: first run everything coarsely and calculate tolerances,
            %then revise
            %TODO: maybe control fidelity with different modes, run nodes
            %wrt specific origins & dimensions; some modes have parameters
            %
            
            t = start;
            
            while t < stop
                for i = 1:length(n.connections)
                    o = n.connections(i).origin.getOutput();
                    n.connections(i).termination.setInput(o);
                end
                      
                endTime = min(t+n.step, stop);
                
                %TODO: parfor here?
                for i = 1:length(n.nodes)
                    n.nodes{i}.run(t, endTime);
                end
                
                for i = 1:length(n.probes)
                    n.probes(i).collect(endTime);
                end
                
                t = t + n.step;
                sprintf('t = %f', t)
            end
        end

        % Performs simulations in both DIRECT_MODE and DEFAULT_MODE 
        % and plots the results together. Data collected by each Probe 
        % are plotted in a separate figure. 
        % 
        % start: simulation start time
        % stop: simulation end time
        % taus (optional): list of time constants for filtering each probe
        %   (default 0.01 for each)
        function runAndPlot(n, start, stop, varargin)
            taus = NemoUtils.getOptionalArg(varargin, 1, .01, 'points', [length(n.probes) 1]); 
            n.reset();
            
            lastFigNum = max([0; findobj('Type', 'figure')]);
            
            n.setSimulationMode(ModeConfigurable.DEFAULT_MODE);
            n.run(start, stop);
            for i = 1:length(n.probes) 
                h = plotProbe(n.probes(i), taus(i), lastFigNum+i, 1, 'b');
                ls = get(h, 'LineStyle');
                if ~isempty(ls) && ~strcmp(ls, 'none') 
                    set(h, 'Color', [.4 .4 .4]);
                end
            end
            n.reset();
            
            n.setSimulationMode(ModeConfigurable.DIRECT_MODE);
            n.run(start, stop);
            for i = 1:length(n.probes) %TODO: no tau for spikes
                h = plotProbe(n.probes(i), taus(i), lastFigNum+i, 0, 'k');
                ls = get(h, 'LineStyle');
                if ~isempty(ls) && ~strcmp(ls, 'none')
                    set(h, 'LineWidth', 2)
                    legend('spiking', 'direct')
                end
            end
%             n.reset();
        end
        
        % origin: An Origin within this Network that is to be available for 
        %       connection outside this Network (within an enclosing, 'direct', 
        %       Network)
        function exposeOrigin(n, origin)
            n.origins{length(n.origins)+1} = origin;
        end
        
        % termination: A Termination within this Network that is to be available for 
        %       connection outside this Network (within an enclosing
        %       Network)
        function exposeTermination(n, termination)
            n.terminations{length(n.terminations)+1} = termination;
        end

        % Resets the Nodes and Probes in this Network
        function reset(n)
            for i = 1:length(n.nodes)
                n.nodes{i}.reset();
            end

            for i = 1:length(n.probes)
                n.probes(i).reset();
            end
        end

        % A convenience method. Additions (as well as deletions) can also 
        % be performed through the public nodes property. 
        % 
        % node: The Node to add 
        function addNode(n, node)
            n.nodes{length(n.nodes)+1} = node;
        end

        % node: Node to remove from Network
        function removeNode(n, node)
            toRemove = [];
            for i = 1:length(n.nodes)
                if eq(node, n.nodes{i})
                    toRemove = [toRemove i];
                end
            end
            toKeep = setdiff(1:length(n.nodes), toRemove);
            n.nodes = n.nodes(toKeep);
        end
        
        % A convenience method. Additions (as well as deletions) can also 
        % be performed through the public connections property. 
        % 
        % origin: Source of the connection
        % termination: Destination of the connection
        function c = addConnection(n, origin, termination)
            c = Connection(origin, termination);
            n.connections = [n.connections; c];
        end
        
        % c: connection to be removed from the network
        function removeConnection(n, c)
            toRemove = [];
            for i = 1:length(n.connection)
                if eq(p, n.connections(i))
                    toRemove = [toRemove i];
                end
            end
            toKeep = setdiff(1:length(n.connections), toRemove);
            n.connections = n.connections(toKeep);
        end
        
        % probeable: A Probeable for which to record state history
        % stateName: Name of the state variable to collect
        function probe = addProbe(n, probeable, stateName)
            probe = Probe(probeable, stateName);
            n.probes = [n.probes; probe];
        end
        
        % Sets simulation modes of all contained Nodes
        function setSimulationMode(n, mode)
            setSimulationMode@ModeConfigurable(n, mode);
            for i = 1:length(n.nodes)
                n.nodes{i}.setSimulationMode(mode);
            end
        end
    end
end
