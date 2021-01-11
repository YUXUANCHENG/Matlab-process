classdef CLI
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        all_mean_cal_A = [];
        order_per = [];
        ifjammed = [];
        mean_p = [];
        var_p = [];
        msd = {};
        ISF_en = {};
        tao_en = {};
        T1_cells_list = {};
        v0 = [];
        clog = [];
        v0_en = {};
        vel_m = {};
        dif = [];
        basefolder;
        periodic = 1;
    end
    
    methods
        function obj = CLI(folder)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.basefolder = folder;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

