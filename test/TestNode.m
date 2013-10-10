% A dummy node that multiplies its input by 2
classdef TestNode < Node
    
    methods (Access = public)
        
        function tn = TestNode(name)
            tn.terminations{1} = Termination(name, .01, 1, 1);
            tn.terminations{1}.node = tn;
            tn.origins{1} = Origin('output', 1);
            tn.origins{1}.node = tn;
        end
        
        function run(tn, start, stop, varargin)
            run(tn.terminations{1}, start, stop);
            val = getOutput(tn.terminations{1})*2;
            setOutput(tn.origins{1}, stop, val);
        end
        
    end
end