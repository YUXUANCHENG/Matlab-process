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
                trial.plotInitial();
                 trial.readMDdata();
                 %trial.setMaxFrames(10000);
                 trial.plotLastFrame(2);
                %trial.showVideo(20);
                %trial.saveVideo(50);    
                trial.createCalculator();
                trial.plotTrajectory(50);
                trial.plotVelDistribution();
                trial.readPhi();
                trial.plotRotationVsTranslaion();
                trial.plotCalADistribution();
                trial.cal_msd();
                trial.plotMSD();
                trial.cal_ISF();
                trial.plotISF();
%                 trial.calculator.cal_c_pos();
%                 trial.verifyISF(3);
%                 trial.cleanUp();
%                  trial.readTao();
%                  trial.plotTaoData();
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
        
        function readTao(obj,mode, phi0, mu, del, eScale)
            energyScale = eScale;
            figure(3); hold on; box on;
            set(gcf,'color','w');
            phi = 0.7:0.015:0.7+9*0.015;
            for t_index_i = 1 :9
                v0_file = obj.basefolder + int2str(t_index_i) + "/" + "v0.txt";
                length_file = obj.basefolder + int2str(t_index_i) + "/" + "length_jammed_" + int2str(t_index_i) + ".txt";
                phi_file = obj.basefolder + int2str(t_index_i) + "/" + "phi_" + int2str(t_index_i) + "0.txt";
                try
                v = csvread(v0_file);
                v = sortrows(v, 8);
                lengthS = csvread(length_file);
                NC = length(lengthS) - 2;
                lengthS = lengthS(end);
                current_phi = phi(t_index_i+1);
                try       
                    p_file = csvread(phi_file);
                    current_phi = mean(p_file,'all');
                end
                area = current_phi* lengthS * lengthS / NC;
                scale = sqrt(area/pi);
                %scale = 1;
                T = (v(:,8)./scale).^2;
                T = T./energyScale;
                tao = v(:,7);
                tao = tao * sqrt(energyScale);
                if mode == 1
                    plot(1./T, tao .* sqrt(T));
                elseif mode == 0
                    scatter(abs(current_phi-phi0)^(2/mu)./(T), abs(current_phi-phi0)^(del)*log10(tao .* sqrt(T)));
                elseif mode == 2
                    plot(1./T, tao);
                end
                catch e
                    disp(e);
                end
            end
            ax = gca;
            ax.XScale = "log";
            ax.YScale = "log";
            xlabel('1/T');ylabel('tao * sqrt(T)');
        end
        
        function Angell(obj,mode)
            figure(3); hold on; box on;
            set(gcf,'color','w');
            phi = [];
            T=[];
            tao=[];
            phi_0 = linspace(0.2,0.8337,40);
            for t_index_i = 0 :39
                v0_file = obj.basefolder + int2str(t_index_i) + "/" + "v0.txt";
                phi_file = obj.basefolder + int2str(t_index_i) + "/" + "phi_" + int2str(t_index_i) + "0.txt";
                length_file = obj.basefolder + int2str(t_index_i) + "/" + "length_jammed_" + int2str(t_index_i) + ".txt";
                try
                v = csvread(v0_file);
                v = sortrows(v, 8);
                length = csvread(length_file);
                length = length(end);
                p_file = csvread(phi_file);
                phi = [phi, mean(p_file,'all')];
                %phi = [phi,phi_0(t_index_i+1)];
                %v = v(1:end-1,:);
                T = [T, (v(:,8)./length).^2];
                tao = [tao, v(:,7)];
                catch
                    disp('no v0');
                end
            end
            threash=10000;
            [~,idx]=min(abs(tao-threash));
            if mode == 0
                phi = phi(1:idx)./phi(idx);
                tao = tao(1:idx);
            end
            plot(phi,tao);
            ax = gca;
            ax.YScale = "log";
            xlabel('phi');ylabel('tao');
            figure(4); hold on; box on;
            set(gcf,'color','w');
            plot(T);
        end
    
        function Arrhenius(obj,mode)
            threash=10000;
            energyScale = 1;
            figure(3); hold on; box on;
            set(gcf,'color','w');
            phi = 0.7:0.015:0.7+9*0.015;
            for t_index_i = 8:9
                v0_file = obj.basefolder + int2str(t_index_i) + "/" + "v0.txt";
                length_file = obj.basefolder + int2str(t_index_i) + "/" + "length_jammed_" + int2str(t_index_i) + ".txt";
                phi_file = obj.basefolder + int2str(t_index_i) + "/" + "phi_" + int2str(t_index_i) + "0.txt";
                try
                v = csvread(v0_file);
                v = sortrows(v, 8, 'descend');
                lengthS = csvread(length_file);
                NC = length(lengthS) - 2;
                lengthS = lengthS(end);
                current_phi = phi(t_index_i+1);
                try       
                    p_file = csvread(phi_file);
                    %current_phi = mean(p_file,'all');
                end
                area = current_phi* lengthS * lengthS / NC;
                scale = sqrt(area/pi);
                %scale = 1;
                T = (v(:,8)./scale).^2;
                T = T./energyScale;
                tao = v(:,7);
                tao = tao * sqrt(energyScale);
                
                [~,idx]=min(abs(tao-threash));
                if mode == 0
                    T = T(1:idx)./T(idx);
                    tao = tao(1:idx);
                end
                plot(1./T,tao);
                ax = gca;
                ax.YScale = "log";
                xlabel('Tg/T');ylabel('tao');
                
                catch e
                    disp(e);
                end
            end
         
        end
        
    end
end

