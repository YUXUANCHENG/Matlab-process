classdef CLI_DPM < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        all_mean_cal_A = [];
        order_per = [];
        ifjammed = [];
        mean_p = [];
        var_p = [];
        sysProperty = {};
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
        function obj = CLI_DPM(folder)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.basefolder = folder;
        end
        
        function trial = createTrial(obj,folder, i, j)
            fileReader = FileReader_back();
            trial = Trial_DPM(i,j,folder,fileReader);
        end

        function pipline(obj, trial)
                %trial.plotInitial();
                trial.readMDdata();
                trial.plotLastFrame(2);
                %trial.showVideo(20);
                %trial.saveVideo(50);
                trial.createCalculator();
                trial.plotVelDistribution();
                trial.readPhi();
                %trial.plotRotationVsTranslaion();
                %trial.plotCalADistribution();
                %trial.cal_msd();
                %trial.plotMSD();
                trial.cal_ISF();
                trial.plotISF();
        end       
        
        function compare(obj, folderList)
            for folder = folderList
                trial = obj.createTrial(folder, 9, 5);
                %obj.createTrial(folder, 1, 3);
                obj.pipline(trial);
            end
        end
        
        function plotScalling(obj)
            figure(3); hold on
            for i = 1:10
                T = [];
                tao = [];
                for j = 1:10
                    T = [T, obj.sysProperty{i,j}.Temp];
                    tao = [tao, obj.sysProperty{i,j}.tao];
                end
                plot(1./T,tao.*sqrt(T));
            end
            ax = gca;
            ax.XScale = "log";
            ax.YScale = "log";
            xlabel('1/v^2');ylabel('tao * v');
        end
        
        function readSysProperty(obj,index_i, index_j, index_i_s, index_j_s)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for t_index_i = index_i_s: index_i
                for t_index_j = index_j_s: index_j
                    try
                        trial = obj.createTrial(obj.basefolder, t_index_i,t_index_j);
                        obj.pipline(trial);
                        obj.sysProperty{t_index_i + 1, t_index_j + 1} = trial;
                    catch e
                        fprintf(1,"%s", e.message);
                        continue
                    end
                end
            end
        end
        
    end
end

