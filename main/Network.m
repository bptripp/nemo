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
    
    properties (SetAccess = private)
        % variables for efficiently running direct mode and population mode
        pmOutputSize = []; %total number of output dimensions from all origins
%         pmTermNodes = []; %list of node indices for each dimension of termination output
        pmPSCInds = []; % indices of PSC state that feed a population's state (state-dim*#terminations x 1)
        pmMap = []; %sparse matrix that maps origins to PSC filter states
        %pmState = []; %list of all state variables in the network
        %pmOutput = []; %list of output of all origin dimensions
        pmDims = []; %list of dimension of state of each node (non-zero for Populations)
%         pmState = [];
        pmTau = [];
        pmFun = []; %list of functions of state associated with each origin 
        pmFunInd = [];
%         pmBiasFun = [];
%         pmNoiseFun = [];
        pmProbeInd = []; % cell array of output indices recorded by each probe (empty for probes that record something other than origin output)
        pmNodes = []; % cell array of nodes to run explicitly in direct/population mode, including sub-networks and custom nodes
        pmTerm = []; % terminations onto nodes other than Populations (values are passed to these manually)
        pmTermInd = []; % cell array of output indices that feed each manual termination 
        pmOriginInd = []; % cell array of output indices that feed each exposed origin
        pmExpTerm = []; % cell array of exposed terminations per node (used only for Populations)
    end
    
    methods (Access = public)
        
        % step: constant step size (s) 
        function n = Network(step)
            n.step = step;
        end

        % start: beginning of simulation time (s)
        % stop: end of simulation time (s)
        function run(n, start, stop)
        	if (n.simulationMode == ModeConfigurable.DIRECT_MODE || n.simulationMode == ModeConfigurable.POPULATION_MODE)
        		runDirect(n, start, stop, 0);
        	else 
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
        end
        
        % Run in DIRECT_MODE. 
        %
        % start: beginning of simulation time (s)
        % stop: end of simulation time (s) 
        % skipSetup: if zero, builds direct or population mode model before
        %   running; otherwise assumes the existing model is current (use
        %   this to start the simulation faster if you know no changes have
        %   been made)
        function runDirect(n, start, stop, skipSetup)
            if ~skipSetup 
                setUpDirectModel(n);
            end
            pscState = zeros(size(n.pmTau));
            popState = zeros(sum(n.pmDims), 1);
        	output = zeros(n.pmOutputSize, 1);
        	
        	% PSC filter parameters
        	pscInputGains = n.step ./ n.pmTau;
        	pscStateGains = 1-pscInputGains;
        	
%             if n.simulationMode == ModeConfigurable.POPULATION_MODE
%                 %TODO: set up models as needed
%             end

			probeTime = (start+n.step):n.step:stop; 
			probeHistorySize = 0;
			for i = 1:length(n.pmProbeInd)
				probeHistorySize = probeHistorySize + length(n.pmProbeInd{i});
			end
			probeHistory = zeros(probeHistorySize, length(probeTime));
			
            numNodes = length(n.nodes);
            t = start;
            for i = 1:length(probeTime)
                % PSC filter 
            	pscInput = n.pmMap * output; % mimics Termination transforms 
            	pscState = pscState .* pscStateGains + pscInput .* pscInputGains; % mimics Termination timeConstants
                
                % sum PSC states of different terminations into population states
                popStateInd = 0;
                popState = zeros(size(popState));
                for j = 1:numNodes 
                    sd = n.pmDims(j);
                    if sd > 0
                        pscInds = n.pmPSCInds{j}; 
                        popState(popStateInd + (1:sd)') = sum(reshape(pscState(pscInds), sd, numel(pscInds)/sd), 2);
                        
                        % add contributions from exposed terminations
                        exposedTerminations = n.pmExpTerm{j};
                        for k = 1:length(exposedTerminations)
                            run(exposedTerminations{k}, probeTime(i)-n.step, probeTime(i));
                            ot = getOutput(exposedTerminations{k});
                            popState(popStateInd + (1:sd)') = popState(popStateInd + (1:sd)') + ot;
                        end
                    end 
                    
                    popStateInd = popStateInd + sd;
                end
                
                % pass values along connections onto non-Populations
                for j = 1:length(n.pmTerm)
                    n.pmTerm{j}.setInput(output(n.pmTermInd{j}));
                end
                
                % run specialized Nodes (not FunctionInputs or Populations)
                for j = 1:length(n.pmNodes)
                    run(n.pmNodes{j}, probeTime(i)-n.step, probeTime(i));
                end
                
                % calculate output as functions of population state
            	outputInd = 0;
                for j = 1:length(n.pmFun)
            		f = n.pmFun{j};
                    if isa(f, 'function_handle') % for DecodedOrigins on Populations and FunctionInputs we run associated functions directly
                        fState = popState(n.pmFunInd{j});
                        if isempty(fState) %a function of time
                            fState = probeTime(i);
                        end
                        o = f(fState);
                    elseif isa(f, 'Origin') % for Origins on other types of Nodes we just run the Node and get the output directly from the Origin
                        o = getOutput(f);
                    else 
                        error('Bug: this should be a function_handle or Origin');
                    end
                    
            		%TODO: add bias and noise from population mode model as needed
                    
            		output(outputInd+(1:length(o))) = o;
            		outputInd = outputInd + length(o);
                end
                
                % set exposed origin values
                for j = 1:length(n.origins)
                    setOutput(n.origins{j}, probeTime(i), output(n.pmOriginInd{j}));
                end
                
                % save probe data 
            	probeInd = 0;
                for j = 1:length(n.pmProbeInd)
                    ind = n.pmProbeInd{j};
                    probeHistory(probeInd+(1:length(ind)),i) = output(ind);
                    probeInd = probeInd + length(ind);
                end
            	
                t = t + n.step;
                if ~mod(i, 100)
                    fprintf('t = %f\n', t)
                end
            end
            
            % set histories (data for origins, zeros for everything else) 
            probeInd = 0;
            for i = 1:length(n.probes)
                ind = n.pmProbeInd{i};
                if ~isempty(ind)
                    history = probeHistory(probeInd+(1:length(ind)),:);
                    setHistory(n.probes(i), probeTime, history);
                    probeInd = probeInd + length(ind);
                end
            end
        end
        
        % Builds a compact model of a network for fast DIRECT_MODE
        % simulations. 
        function setUpDirectModel(n)
            n.pmOutputSize = 0;
            n.pmPSCInds = cell(size(n.nodes)); %cell array of PSC state indices for each node 
            n.pmMap = [];
            n.pmDims = zeros(length(n.nodes), 1);
            n.pmTau = [];
            n.pmFun = cell(1, 0); 
            n.pmFunInd = cell(1,0);
            n.pmTerm = cell(1,0);
            n.pmTermInd = cell(1,0);
            n.pmExpTerm = cell(1,0);
            
            originUIDs = cell(1,0); %list of node-name:origin-name in order of input to pmMap
            termUIDs = cell(1,0); %likewise for terminations

            for i = 1:length(n.nodes)
                origins = n.nodes{i}.origins;
                for j = 1:length(origins)
                    UID = getOriginUID(origins{j});
                    for k = 1:origins{j}.dim
                        originUIDs{length(originUIDs)+1} = UID;
                    end
                end
                
                n.pmExpTerm{i} = cell(1,0);
                
                if isa(n.nodes{i}, 'Population')
                    % build list of state size for each node
                    n.pmDims(i) = length(n.nodes{i}.radii);

                    for j = 1:length(origins)
                        if isa(origins{j}, 'DecodedOrigin')
                            ind = length(n.pmFun)+1;
                            n.pmFun{ind} = origins{j}.f;
                            n.pmFunInd{ind} = sum(n.pmDims(1:i-1)) + (1:n.pmDims(i));  % indices of pop state vector that make up input to this function
                        end
                    end    
                    
                    % we only aggregate Terminations on Populations 
                    terminations = n.nodes{i}.terminations;
                    n.pmPSCInds{i} = [];
                    for j = 1:length(terminations)
                        UID = getTerminationUID(terminations{j});
                        ind = length(termUIDs) + (1:size(terminations{j}.transform, 1))';
                        n.pmPSCInds{i} = [n.pmPSCInds{i} ind]; 

                        % build UID list with each input dim
                        for k = 1:size(terminations{j}.transform, 1)
                            termUIDs{length(termUIDs)+1} = UID;
                        end
                        
                        for k = 1:length(n.terminations)
                            exposedUID = getTerminationUID(n.terminations{k});
                            if strcmp(UID, exposedUID)
                                n.pmExpTerm{i}{length(n.pmExpTerm{i})+1} = n.terminations{k};
                            end
                        end
                    end
                    
                elseif isa(n.nodes{i}, 'FunctionInput')
                    ind = length(n.pmFun)+1;
                    n.pmFun{ind} = n.nodes{i}.func;
                    n.pmFunInd{ind} = []; % empty list of input states indicates function of time
                else 
                    n.pmNodes{length(n.pmNodes)+1} = n.nodes{i};
                    for j = 1:length(origins)
                        ind = length(n.pmFun)+1;
                        n.pmFun{ind} = origins{j}; % we will get values straight from the origin in this case
                    end
                end

            end

            n.pmOutputSize = length(originUIDs);
            n.pmMap = zeros(length(termUIDs), length(originUIDs));
            n.pmTau = ones(length(termUIDs), 1); %avoid divide by 0 if not set
            for i = 1:length(n.connections)
                termination = n.connections(i).termination;
                
                originUID = getOriginUID(n.connections(i).origin);
                originInd = find(strcmp(originUIDs, originUID)); 
                
                if isa(termination.node, 'Population') 
                    % build pmMap ...
                    termUID = getTerminationUID(termination);
                    termInd = find(strcmp(termUIDs, termUID));
                    n.pmMap(termInd, originInd) = termination.transform;

                    % collect time constants ...
                    n.pmTau(termInd) = termination.timeConstant;
                else % connections that terminate onto other node types are run directly
                    ind = length(n.pmTerm)+1;
                    n.pmTerm{ind} = termination;
                    n.pmTermInd{length(n.pmTerm)} = originInd;
                end

            end
            
            % setup for passing exposed origin values out
            for i = 1:length(n.origins)
                originUID = getOriginUID(n.origins{i});
                originInd = find(strcmp(originUIDs, originUID)); 
                n.pmOriginInd{i} = originInd;
            end
            
            % setup for getting exposed termination values for Populations
            % (custom nodes are run normally so they don't have to be
            % handled here)
            for i = 1:length(n.terminations)
                termUID = getTerminationUID(n.terminations{i});
                termInd = find(strcmp(termUIDs, termUID));
            end

            n.pmProbeInd = cell(1, length(n.probes));
            for i = 1:length(n.probes)
                p = n.probes(i).probeable;
                if isa(p, 'Origin')
                    if ~isempty(p.node)
                        UID = getOriginUID(p);
                        n.pmProbeInd{i} = find(strcmp(originUIDs, UID));
                    end
                else 
                    n.pmProbeInd{i} = [];
                end
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
            origin.node = n;
        end
        
        % termination: A Termination within this Network that is to be available for 
        %       connection outside this Network (within an enclosing
        %       Network)
        function exposeTermination(n, termination)
            n.terminations{length(n.terminations)+1} = termination;
            termination.node = n;
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

function result = getOriginUID(origin)
    result = [];
    if ~isempty(origin.node) 
        result = cat(2, origin.node.name, ':', origin.name);
    else 
        warning('Origin %s has unknown Node and will be ignored in DIRECT_MODE.', origin.name)
    end
end

function result = getTerminationUID(termination) 
    result = [];
    if ~isempty(termination.node) 
        result = cat(2, termination.node.name, ':', termination.name);
    end
end
