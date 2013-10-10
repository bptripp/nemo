% Termination of a projection that can be simulated in a reduced form. 

% TODO: there should be a common superclass of this and Termination,
%      referenced by Node
% TODO: make probeable

classdef ReducibleTermination < handle
    
    properties (Access = public)
        % A Name that distinguishes this from other Termination on the same
        %   Node
        name = [];
    end
    
    properties (SetAccess = private)
        weights = [];
        U = [];
        S = [];
        V = [];
        
        timeConstant = [];
        stepRatio = 1;
        reducedDim = []; % empty means not reduced
    end
    
    properties (Access = private)
        A = [];
        B = [];
        input = [];
        state = [];
        time = [];
    end
    
    methods (Access = public)
        
        function rt = ReducibleTermination(name, timeConstant, weights, stepRatio)
            rt.name = name;
            rt.timeConstant = timeConstant;
            rt.stepRatio = stepRatio;

            rt.weights = weights;            
            [rt.U, rt.S, rt.V] = svd(weights);
        end
        
        function setReducedDim(rt)
            %TODO
        end
        
        % Called by the Network at the beginning of each time step
        % 
        % input: The output of a connected Origin (e.g. unit spike
        %   impulses; decoded variables)
        function setInput(t, input)   
            assert(size(input,1) == size(t.B, 2), 'Expected an input vector of length %i, got %i', size(t.B, 2), size(input,1));
            assert(size(input,2) == 1, 'Expected a column vector as input');
            t.input = input;
        end

        % output: The result of the Termination dynamics, which enters the
        %   associated Node 
        function output = getOutput(t)
            output = t.state;
        end
        
        % Runs (integrates) the Termination for a certain amount of
        % simulation time, and updates the output. 
        % 
        % time: simulation time for which to integrate the Termination
        function run(t, start, stop)
            time = stop - start;
            subStepTime = time / t.stepRatio;
                        
            ff = t.B * t.input; %feedforward component constant within run
            for i = 1:t.stepRatio
                dsdt = t.A * t.state + ff;
                t.state = t.state + dsdt * subStepTime;
%                 sprintf('input %f, state %f, fb %f, ff %f, dsdt %f', t.input, t.state, t.A*t.state, ff, dsdt)
            end            
            t.time = stop;
        end
        
        % Resets internal state to zero
        function reset(t)
            t.time = [];
            t.state = zeros(size(t.B, 1), 1);
        end
        
    end    
    
end