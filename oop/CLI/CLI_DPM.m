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
                %trial.cal_ISF();
                %trial.plotISF();
                trial.calculator.cal_c_pos();
                trial.verifyISF(3);
                trial.cleanUp();
        end       
        
        function compare(obj, folderList)
            for folder = folderList
                trial = obj.createTrial(folder, 0, 5);
                obj.pipline(trial);
            end
        end
        
        function plotScalling(obj)
            figure(3); hold on
            for i = 1:1
                T = [];
                tao = [];
                for j = 2:10
                    T = [T, obj.sysProperty{i,j}.Temp * (obj.sysProperty{i,j}.lengthscale(end))^2];
                    tao = [tao, obj.sysProperty{i,j}.tao *10^3*100000/5000];
                end
                plot(1./T,tao.*sqrt(T));
            end
            ax = gca;
            ax.XScale = "log";
            ax.YScale = "log";
            xlabel('1/T');ylabel('tao * sqrt(T)');
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
        
        function readTao(obj)
            figure(3); hold on; box on;
            set(gcf,'color','w');
            phi = 0.7:0.02:0.7+9*0.02;
            phi0 = 0.849;
            mu = 1.30;
            %del = 1.95;
            del = 2.9;
            for t_index_i = 0 :9
                v0_file = obj.basefolder + int2str(t_index_i) + "/" + "v0.txt";
                try
                v = csvread(v0_file);
                v = sortrows(v, 8);
                v = v(1:end-1,:);
                T = v(:,8).^2;
                tao = v(:,7);
                plot(1./T, tao .* sqrt(T));
                %scatter(abs(phi(t_index_i)-phi0)^(2/mu)./T, abs(phi(t_index_i)-phi0)^(del)*log(tao .* sqrt(T)));
                catch
                    disp('no v0');
                end
            end
            ax = gca;
            ax.XScale = "log";
            ax.YScale = "log";
            xlabel('1/T');ylabel('tao * sqrt(T)');
        end
        
    end
end

