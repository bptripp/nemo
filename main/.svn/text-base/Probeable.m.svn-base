% An object that has some internal state that can be probed. (See also 
% Probe.) This provides a way of collecting and retaining selected 
% simulation data for analysis after a run. 
classdef Probeable < handle 
    
    methods (Abstract) 

        % names: Cell array of the names of all the state variables of this
        %       object that can be probed 
        names = getStateNames(p);

        % name: Name of a state variable 
        % dim: Dimension of named state variable
        dim = getDimension(p, name);

        % name: Name of a state variable 
        % time: End of most recent time step
        % state: Value of state variable at end of most recent time step
        [time, state] = getState(p, name);        
    end
end