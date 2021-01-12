classdef FileReader_back < FileReader
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       
    end
    
    methods
        function obj = setFolder(obj, basefolder, t_index_i, t_index_j)
            obj.folder = basefolder + int2str(t_index_i) + "/";
            obj.t_index_i = t_index_i;
            obj.t_index_j = t_index_j;
        end
    end
end

