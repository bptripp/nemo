% Source of output of a Node (e.g. the axon of a cell, including some number 
% of presynaptic specializations, or a group of axons that correspond to a 
% population, or a linear decoding of population output).
classdef Origin < Probeable & ModeConfigurable
    
    properties (Access = public)
        % A name that distinguishes this origin from others of the same
        % Node, e.g. different decoded outputs. 
        name = [];
        
        % Length of output vector
        dim = [];
    end
    
    properties (Access = private)
        output = [];
    end
    
    properties (Access = private)
        lastTime = [];
    end
        
    methods (Access = public)
        
        % name: A name that distinguishes this origin from others of the same
        %       Node, e.g. different decoded outputs. 
        function o = Origin(name, dim)
            o.name = name;     
            o.dim = dim;
            o.reset();
        end
        
        % output: An o.dim by 1 vector containing the instantaneous value 
        %       of the output of this Origin, normally set in setOutput(...) 
        %       during a run(...) step of the associated Node
        function output = getOutput(o)
            output = o.output;
        end
        
        % This should be called by the associated Node before completion of
        % a run(...), as a way of making Node output available to other
        % Nodes. 
        %         
        % time: Simulation time (at end of run step)
        % output: Value of the output (o.dim by 1)
        function setOutput(o, time, output)
            assert(size(output,1) == o.dim, 'Output has wrong dimension') 
            assert(size(output,2) == 1, 'Output should be a column vector') 
            
            o.output = output;
            o.lastTime = time;
        end
        
        function reset(o)
            o.lastTime = 0;
            o.output = zeros(o.dim, 1);
        end

        % see Probeable
        function names = getStateNames(o) 
            names = {'output'};
        end

        % see Probeable
        function dim = getDimension(o, name)
            assert(strcmp(name, 'output'), 'Unknown state name');
            dim = o.dim;
        end

        % see Probeable
        function [time, state] = getState(o, name)
            assert(strcmp(name, 'output'), 'Unknown state name');
            time = o.lastTime;
            state = o.output;
        end
        
    end
end
