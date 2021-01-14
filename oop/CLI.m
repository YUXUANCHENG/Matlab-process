classdef CLI < handle
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
        v0_en = {};
        vel_m = {};
        dif = [];
        basefolder;
    end
    
    methods
        function obj = CLI(folder)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.basefolder = folder;
        end
    end
    
    methods(Static)     
        function compare_DPM_SP(folderList)
            for folder = folderList
                fileReader = FileReader_back();
                trial = Trial(3,2,folder,fileReader);
                %trial.plotInitial();
                trial.readMDdata();
                trial.plotLastFrame(2);
                trial.createCalculator();
                trial.plotVelDistribution();
                trial.cal_msd();
                trial.plotMSD();
                delete(fileReader);
                delete(trial);
            end
        end
        
    end
end

