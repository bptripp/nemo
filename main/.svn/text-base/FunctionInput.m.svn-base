% An abstract model input, used to provide signals that mimic neurons or 
% sense organs outside the scope of the network. These signals are
% typically explicit functions of time. 
classdef FunctionInput < Node
    
    properties (Access = public)
        func = []; 
    end
    
    methods (Access = public)
        
        % func: n-dimensional function of time
        function in = FunctionInput(func)
            in.func = func;
            class(func)
            assert(strcmp(class(func), 'function_handle'), 'Expected a function handle');
            testOutput = func(0);
            assert(size(testOutput, 2) == 1, 'Function output should be a column vector')
            dim = size(testOutput, 1);
            in.origins{1} = Origin('X', dim);
            in.run(0, 0);
        end

        % see Node
        function run(in, start, stop)
            in.origins{1}.setOutput(stop, double(in.func(stop)));
        end

        % origin: this Node's single Origin, which outputs the function
        %   given in the constructor
        function origin = getOrigin(in)
            origin = in.origins{1};
        end
        
    end
end