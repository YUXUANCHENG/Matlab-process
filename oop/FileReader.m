classdef FileReader < handle
    %FILEREADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v0
        folder
        lengthscale
        coordinate
        t_index_i
        t_index_j
        vel
        cal_A
        phi
        isfData
        endpoint = 35;
        interval = 15;
        flow;
    end
    
    methods
        function setInterval(obj,endpoint,interval)
            obj.endpoint = endpoint;
            obj.interval = interval;
        end
        
        function obj = setFolder(obj, basefolder, t_index_i, t_index_j)
            obj.folder = basefolder + int2str(t_index_i) + "_" + int2str(t_index_j) + "/";
            obj.t_index_i = t_index_i;
            obj.t_index_j = t_index_j;
        end
        
        function readV0(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            v0_file = obj.folder + "v0.txt";
            try
                obj.v0 = csvread(v0_file);
            catch
                disp('no v0');
            end
        end
        
        function [flow1, flow_rate] = readFlowR(obj)
            obj.flow = csvread(obj.folder + "flowRate.txt");
            flow1 = obj.flow;
            limit = obj.endpoint;
%             limit = 19;
%             limit = 38;
%              limit = 59;
            if size(obj.flow,1) > limit
%                 flow = flow(limit + 1:end);
                 flow = obj.flow(limit - obj.interval:limit + 1);
%                 flow = flow(limit - 45:limit + 1);
%                 flow = flow(limit - 20:limit + 1);
%                  flow = flow(limit - 10:limit + 1);
%                 flow = flow(limit+1:end);
                flow_rate = 10*mean(flow,'all');
                   
            else
                flow = obj.flow(end);
                flow_rate = -1;
            end
            plotF = 0;
            if (plotF)
            figure(5); hold on;
            set(gcf,'color','w');
%             edges = linspace(0, 0.04, 40);
%             h1=histogram(10*flow, 'BinEdges', edges);
%             h1.BinWidth = 0.005;
%             h1.BinCount = 10;
            interval = csvread(obj.folder + "flowInterval.txt");
%             interval = interval(floor(size(interval,1)/2):end);
            interval = interval(30:end);
            edges = linspace(0, 9000, 18);
%              [~,edges] = histcounts(log10(interval));
            h1=histogram(interval, 'BinEdges', edges, 'Normalization','probability');
            height = h1.Values;
%             h1=histogram(interval, 'BinEdges', 10.^edges, 'Normalization','probability');
            set(gca,'YScale','log');           
%              set(gca, 'xscale','log');
            figure(6); hold on; box on;
            set(gcf,'color','w');
            plot(edges(1:end-1),height)
            set(gca,'YScale','log');  
            set(gca,'FontSize',15)
            xlabel('$\delta t$','Interpreter','latex');
            ylabel("$P$",'Interpreter','latex');
            end
        end
        
        function friction = readFriction(obj)
            friction = csvread(obj.folder + "flowRate.txt");
        end
        
        function readInitial(obj)
            extend = "_jammed_" + int2str(obj.t_index_i) +".txt";
            coordinate_file = obj.folder + "jam" + extend;
            length_file = obj.folder + "length" + extend;
            cal_A_file = obj.folder + "calA" + extend;
            contact_file = obj.folder + "contact" + extend;
            
            try
                packing_contact = dlmread(contact_file);
                packing_contact = packing_contact(end,:);
            catch
            end
            obj.coordinate = csvread(coordinate_file);
            %length_file = folder + "length" + extend;
            obj.lengthscale = csvread(length_file);
            obj.cal_A = csvread(cal_A_file);
        end
        
        function readMDdata(obj)
            if isempty(obj.lengthscale)
                obj.readInitial()
            end
            extend1 = "_" + int2str(obj.t_index_i) + int2str(obj.t_index_j) +".txt";
            coordinate_file = obj.folder + "jam" + extend1;
            obj.coordinate = csvread(coordinate_file);
            phi_file = obj.folder + "phi" + extend1;
            obj.phi = csvread(phi_file);
            try
                cal_A_file = obj.folder + "calA" + extend1;
                obj.cal_A = csvread(cal_A_file);
            catch
                disp('no calA file')
            end
%             try
%                 contact_file = obj.folder + "contact" + extend1;
%                 contact = dlmread(contact_file);
%             catch
%                 disp('no contact file')
%             end
            v_file = obj.folder + "v" + extend1;
            obj.vel = csvread(v_file);
        end
        
        function readTao(obj)
            extend = "_" + int2str(obj.t_index_i) + int2str(obj.t_index_j) +".txt";
            ISF_file = obj.folder + "isf" + extend;
            obj.isfData = csvread(ISF_file);
        end
    end
end

