% A point of input to a Node (e.g. a post-synaptic specialization). The
% default implementation passed input values through first-order
% exponential dynamics, which serve as a simple model of post-synaptic 
% current dynamics. Subclass Termination to change this. 
classdef Termination < Probeable
    
    properties (Access = public)
        % A Name that distinguishes this from other Termination on the same
        %   Node
        name = [];
    end
    
    properties (SetAccess = private)
        transform = [];
        timeConstant = [];
        stepRatio = 1;
        bias = [];
        biasEncoders = [];
    end
    
    properties (Access = private)
        A = [];
        B = [];
        input = [];
        state = [];
        time = [];
    end
    
    methods (Access = public)

        % name: A Name that distinguishes this from other Termination on 
        %   the same Node        
        % timeConstant: time constant of exponential decay
        % transform: a matrix that maps input to state variables 
        % stepRatio: number if internal steps per network step (>=1)
        function t = Termination(name, timeConstant, transform, stepRatio)
            t.name = name;
            t.transform = transform;
            t.timeConstant = timeConstant;
            t.stepRatio = stepRatio;
            t.bias = zeros(size(transform, 1), 1);

            assert(stepRatio >= 1, 'stepRatio must be >=1')
            t.stepRatio = round(stepRatio);
            
            t.A = -eye(size(transform,1)) / timeConstant;
            t.B = transform ./ timeConstant;
            reset(t);
        end

        % bias: Static bias in the output of this Termination, which modifies
        %   baseline neuron activity. It is convenient to include a bias term
        %   here because in some cases bias should be added / removed from a
        %   model in conjunction with a Termination (particularly when
        %   enabling / disabling a Parisien projection). 
        function setBias(t, bias)
            assert(size(bias, 1) == size(t.transform, 1), 'Expected a column vector same size as transform output')
            assert(size(bias, 2) == 1, 'Expected a column vector')
            t.bias = bias;
        end

        % biasEncoders: Encoders for each neuron in the associatd population
        %   that are specific to this Termination (allows encoding of
        %   Parisien bias signals orthogonally to represented signals)
        function setBiasEncoders(t, biasEncoders)
            t.biasEncoders = biasEncoders;
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
            output = t.state + t.bias;
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

        % see Probeable
        function names = getStateNames(t) 
            names = {'state'};
        end

        % see Probeable
        function dim = getDimension(t, name)
            assert(strcmp(name, 'state'), 'No such state')
            dim = size(t.A,1);
        end

        % see Probeable
        function [time, state] = getState(t, name)
            assert(strcmp(name, 'state'), 'No such state')
            time = t.time;
            state = t.state;
        end
    end
end