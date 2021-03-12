classdef CLI_Disk < CLI_DPM
    %CLI_DISK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods
        function createTrial(obj,folder, i, j)
            obj.fileReader = FileReader_back();
            obj.trial = Trial_Disk(i,j,folder,obj.fileReader);
        end
    end
end

