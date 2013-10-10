classdef NemoUtils < handle
    
    methods (Static)
        
        % Extracts elements from a varargin cell array.         
        % 
        % varargin: cell array of optional arguments to some other method
        % index: the index of the optional argument of interest
        % default: a default value if the argument does not exist (is empty
        %   or the list is too short)
        % name: name of the argument (for error messages)
        % dims: expected size of argument (row vector; use [] if unknown). 
        %   If the argument of interest has one element, it is expanded to 
        %   a uniform matrix of this size
        function result = getOptionalArg(varargin, index, default, name, dims)
            assert(isempty(dims) || size(dims, 1) == 1, 'dims should be a row vector')
            
            if length(varargin) >= index && ~isempty(varargin{index})
                result = varargin{index};

                if numel(result) > 1
                    for i = 1:length(dims)
                        assert(size(result, i) == dims(i), sprintf('Size of %s should be %i in dimension %i', name, dims(i), i));
                    end
                end
            else 
                result = default;
            end
            
            if numel(result) == 1 && ~isempty(dims)
                result = result * ones(dims);
            end
        end

    end
end