% Something that has runs in different simulation modes. Simulations modes
% trade off between performance and realism. 
classdef ModeConfigurable < handle
    
    properties (Constant = true)
        % A low-cost mode in which it is not necessary to simulate neurons.
        % Populations are taken to represent variables exactly and origins
        % are taken to decode functions exactly. 
        DIRECT_MODE = 1; 
        
        % Neurons do not spike, but communicate in terms of firing rates.
        % Firing rates have no dynamics, i.e. the amount of time for which 
        % a neuron runs is not considered. 
        % Variables are decoded from firing rates. 
        CONSTANT_RATE_MODE = 2;
        
        % Neurons do not spike, but communicate in terms of firing rates. 
        % Variables are decoded from firing rates. 
        RATE_MODE = 3;
        
        % Neurons communicate through spikes, and represented variables are
        % decoded from spiking activity.         
        DEFAULT_MODE = 4;
        
        % Principal components of neuron responses are simulated. 
        PC_MODE = 5;
    end
    
    properties (SetAccess = private)
        simulationMode  = ModeConfigurable.DEFAULT_MODE;        
    end
    
    properties (Access = public)
        nPC = 10;
    end
    
    methods (Access = public)
        
        % mode: One of the constant DIRECT_MODE, RATE_MODE, or DEFAULT_MODE
        function setSimulationMode(mc, mode)
            assert(mode == ModeConfigurable.DIRECT_MODE ...
                || mode == ModeConfigurable.RATE_MODE ...
                || mode == ModeConfigurable.DEFAULT_MODE ...
                || mode == ModeConfigurable.CONSTANT_RATE_MODE, 'Unknown simulation mode')
            mc.simulationMode = mode;
        end
    end
    
end