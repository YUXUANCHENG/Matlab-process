classdef CLI_DPM < handle
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
        fileReader;
        trial;
    end
    
    methods
        function obj = CLI_DPM(folder)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.basefolder = folder;
        end
        
        function createTrial(obj,folder, i, j)
            obj.fileReader = FileReader_back();
            obj.trial = Trial_DPM(i,j,folder,obj.fileReader);
        end

        function pipline(obj)
                obj.trial.plotInitial();
                obj.trial.readMDdata();
                obj.trial.plotLastFrame(2);
                %obj.trial.showVideo(20);
                %obj.trial.saveVideo(50);
                obj.trial.createCalculator();
                obj.trial.plotVelDistribution();
                obj.trial.plotRotationVsTranslaion();
                obj.trial.plotCalADistribution();
                obj.trial.cal_msd();
                obj.trial.plotMSD();
                delete(obj.fileReader);
                delete(obj.trial);
        end
        
        function compare(obj, folderList)
            for folder = folderList
                obj.createTrial(folder, 1, 0);
                %obj.createTrial(folder, 1, 3);
                obj.pipline();
            end
        end
        
    end
end

