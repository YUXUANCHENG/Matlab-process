classdef CLI_hopper < CLI_DPM
    %CLI_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hopperProperty = []
        a
        b
    end
    
    methods
        function calHopperProperty(obj,index_i, index_j, mode, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.HopperProperty(index_i, index_j, 0, 0, mode, varargin{:})
        end
        
        function HopperProperty(obj,index_i, index_j, index_i_s, index_j_s, mode, varargin)
            %METHOD1 Summary of this method goes here
            obj.hopperProperty = [];
            %   Detailed explanation goes here
            for t_index_i = index_i_s: index_i
                for t_index_j = index_j_s: index_j
                    try
                        fileReader = FileReader();
                        trial = Trial_hopper(t_index_i,t_index_j,obj.basefolder,fileReader);
                        %trial.plotInitial();
                        trial.readV0();
                        trial.readMDdata();
                        trial.plotLastFrame(2);
                        trial.showVideo(50);
                        %trial.saveVideo(50);
                        trial.createCalculator();
                        trial.printCellCount();
                        trial.plotEk();
                        trial.flowRate(mode,varargin);
                        trial.calDensity(20,3);
                        trial.plotDensity();
                        obj.hopperProperty = [obj.hopperProperty, trial];
                    catch e
                        fprintf(1,"%s", e.message);
                        delete(fileReader);
                        continue
                    end
                    trial.cleanUp();
                end
            end
        end
        
        function plotFlowRate(obj)
            flowRate = [];
            error = [];
            width = [];
            for i = 1: length(obj.hopperProperty)
                width = [width, obj.hopperProperty(i).width];
            end
            width = unique(width);
            for w = width
                rate_temp = [];
                for i = 1: length(obj.hopperProperty)
                    if obj.hopperProperty(i).width == w
                        rate_temp = [rate_temp, obj.hopperProperty(i).flow_rate];
                    end
                end
                flowRate = [flowRate, mean(rate_temp)];
                error = [error, std(rate_temp)/sqrt(length(rate_temp))];
            end
            figure(3);
            errorbar(width,flowRate,error,'o')
            %errorbar(1 - gam/(clog(1,1)*lengthscale(1,1)*coordinate(1,3)),clog_p,error,'o')
            xlabel("width");
            ylabel("flow rate");
        end
%         
%         function plotFlowRateVSN_W(obj)
%             flowRate = [];
%             N_W = [];
%             for i = 1: length(obj.hopperProperty)
%                 flowRate = [flowRate, obj.hopperProperty(i).flow_rate{1,1}];
%                 N_W = [N_W, obj.hopperProperty(i).flow_rate{1,2}/obj.hopperProperty(i).width];
%             end             
%             figure(3);
%             scatter(N_W,flowRate,'o')
%             xlabel("N/W");
%             ylabel("flow rate");              
%         end
%         
        function plotFlowRateVSN_W(obj)
            width = [];
            for i = 1: length(obj.hopperProperty)
                width = [width, obj.hopperProperty(i).width];
            end
            width = unique(width);
            %for w = width
            for w = width
                flowRate = [];
                error = [];
                N_left = [];
                for i = 1: length(obj.hopperProperty)
                    if obj.hopperProperty(i).width == w
                        N_left = [N_left, obj.hopperProperty(i).flow_rate{1,2}];
                    end
                end
                N_left = unique(N_left);
                N_left = N_left(N_left > 0);
                for number = N_left
                    rate_temp = [];
                    for i = 1: length(obj.hopperProperty)
                        if obj.hopperProperty(i).width == w
                            search = find(obj.hopperProperty(i).flow_rate{1,2} == number);
                            if ~isempty(search)
                                rate_temp = [rate_temp, obj.hopperProperty(i).flow_rate{1,1}(search)];
                            end
                        end
                    end
                    flowRate = [flowRate, mean(rate_temp)];
                    error = [error, std(rate_temp)/sqrt(length(rate_temp))];
                end   
                figure(3);hold on;
                errorbar(N_left/w,flowRate/w,error,'o')
                %scatter(N_left/w,flowRate/w,'o')
                %plot(N_left,flowRate,'o')
                %errorbar(1 - gam/(clog(1,1)*lengthscale(1,1)*coordinate(1,3)),clog_p,error,'o')
                xlabel("N/width");
                ylabel("flux");  
            end
            legend(string(width));
        end
        
        function plotFlowRateWithN(obj)
            width = [];
            for i = 1: length(obj.hopperProperty)
                width = [width, obj.hopperProperty(i).width];
            end
            width = unique(width);
            for w = width
                flowRate = [];
                error = [];
            %for w = width(end)
                N_left = [];
                for i = 1: length(obj.hopperProperty)
                    if obj.hopperProperty(i).width == w
                        N_left = [N_left, obj.hopperProperty(i).flow_rate{1,2}];
                    end
                end
                N_left = unique(N_left);
                for number = N_left
                    rate_temp = [];
                    for i = 1: length(obj.hopperProperty)
                        if obj.hopperProperty(i).width == w
                            search = find(obj.hopperProperty(i).flow_rate{1,2} == number);
                            if ~isempty(search)
                                rate_temp = [rate_temp, obj.hopperProperty(i).flow_rate{1,1}(search)];
                            end
                        end
                    end
                    flowRate = [flowRate, mean(rate_temp)];
                    error = [error, std(rate_temp)/sqrt(length(rate_temp))];
                end                
            
                figure(3);hold on;
                errorbar(N_left,flowRate,error,'o')
                %errorbar(1 - gam/(clog(1,1)*lengthscale(1,1)*coordinate(1,3)),clog_p,error,'o')
                xlabel("N left inside");
                ylabel("flow rate");   
            end
            legend(string(width));
        end
        
        function plotClogP(obj)
            clog_p = [];
            error = [];
            width = [];
            for i = 1: length(obj.hopperProperty)
                width = [width, obj.hopperProperty(i).width];
            end
            width = unique(width);
            for w = width
                results = [];
                for i = 1: length(obj.hopperProperty)
                    if obj.hopperProperty(i).width == w
                        results = [results, obj.hopperProperty(i).result];
                    end
                end
                clog_p = [clog_p, mean(results)];
                error = [error, std(results)/sqrt(size(results, 2))];
            end
            fitfun = fittype( @(a, b, x) 1 ./ (1 + exp( (x - a) ./ b)));
            %fitted = fit( (obj.logindex(half:end))', abs(obj.ISF(obj.logindex(half:end))), fitfun, 'StartPoint', [10 * (10 - obj.t_index_j),1]);
            fitted = fit(width', clog_p', fitfun, 'StartPoint', [1,0.1]);
            disp({'a = ', fitted.a, ' b = ', fitted.b})
            obj.a = fitted.a;
            obj.b = fitted.b;
            %obj.tao = abs(fitted.tao);
            figure(4); hold on; box on;
            set(gcf,'color','w');
            errorbar(width,clog_p,error,'o')
            %errorbar(1 - gam/(clog(1,1)*lengthscale(1,1)*coordinate(1,3)),clog_p,error,'o')
            plot(fitted);
            xlabel("width");
            %xlabel("alpha");
            ylabel("clogging probability");
        end
        
        function compare(obj, folderList)
            a_list = [];
            b_list = [];
            for folder = folderList
                obj.basefolder = folder;
                obj.calHopperProperty(99, 39, 4, 50);
                obj.plotClogP();
                a_list = [a_list, obj.a];
                b_list = [b_list, obj.b];
            end
            figure(5);
            scatter([10,1,0.1,0.01],a_list,25,'filled');
            xlabel('kb');ylabel('a');
            ax = gca;
            ax.XScale = "log";
            ax.YScale = "Linear";
        end
    end
end

