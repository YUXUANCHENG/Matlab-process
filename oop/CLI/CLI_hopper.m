classdef CLI_hopper < CLI_DPM
    %CLI_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hopperProperty = {}
        a
        b
        angles = [];
        start = 0.5;
        interval = 0.1;
        fInter = 15;
        fEnd = 35;
    end
    
    methods
        
        function setInterval(obj,start,interval,fEnd,fInter)
            obj.start = start;
            obj.interval = interval;
            obj.fEnd = fEnd;
            obj.fInter = fInter;
        end
        
        function calHopperProperty(obj,index_i, index_j, mode, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.HopperProperty(index_i, index_j, 0, 0, mode, varargin{:})
        end
        
        function HopperProperty(obj,index_i, index_j, index_i_s, index_j_s, mode, varargin)
            %METHOD1 Summary of this method goes here
            obj.hopperProperty = [];
            obj.angles = [];
            %   Detailed explanation goes here
            for t_index_i = index_i_s: index_i
                for t_index_j = index_j_s: index_j
                    try
%                         if mod(t_index_j,10)>0
%                             continue;
%                         end
                        fileReader = FileReader();
                        fileReader.setInterval(obj.fEnd,obj.fInter);
                        trial = Trial_hopper(t_index_i,t_index_j,obj.basefolder,fileReader);
%                         trial.plotInitial();
                        trial.width = roundn((obj.start + obj.interval * t_index_j),-4);
%                          trial.width = roundn((0.5 + 1 * t_index_j),-4);
%                            trial.width = roundn((3 + 1 * t_index_j),-4);
%                             trial.width = roundn((3 + 0.4 * t_index_j),-4);
%                         trial.width = roundn((0.5 + 0.1 * t_index_j)*1.04,-4);
%                         trial.width = roundn(3 + 0.1 * t_index_j,-4);
%                         trial.width = roundn(0.3 + 0.06 * t_index_j,-4);
%                         trial.readV0();
%                         trial.readAndPlotFriction();
                           %trial.readMDdata();
                           trial.readFlowR();
                           if (trial.flow_rate < 0)
                               continue
                           end
%                            trial.plotFlow();
%                         trial.cal_eff_width_scale();
%                         if (trial.result == 1) && (trial.width>1.5)
% %                        if (trial.result == 0)&& (trial.width < 1)
% % %                         if (trial.width>8)
%                          trial.plotInitial();
%                          trial.createCalculator();
%                          trial.readMDdata();
% %                          trial.volumnRate();
%                          angle = trial.calculator.cal_vel_distri_in_space(1);
%                          obj.angles = [obj.angles, angle];
%                          trial.plotLastFrame(2);
% % % % % % % % % % %                          trial.createCalculator();
% % % % % % % % %                          trial.printCellCount();
% %                          trial.saveVideo(50);
% % % % % %                         end
%                        trial.showVideo(50);
% % %                         trial.saveVideo(50);
%                          trial.createCalculator();
% % % %                         trial.printCellCount();
%                         trial.plotEk();
%                         trial.flowRate(mode,varargin);
%                         trial.calDensity(20,3);
%                         trial.plotDensity();
                        trial.cleanUp();
%                         obj.hopperProperty = [obj.hopperProperty, trial];
                        obj.hopperProperty{end+1} = trial;
                    catch e
%                         trial.result = 1;
%                         obj.hopperProperty{end+1} = trial;
                        fprintf(1,"%s", e.message);
                        delete(fileReader);
                        continue
                    end
                    trial.cleanUp();
                end
            end
            obj.hopperProperty = [obj.hopperProperty{:}];
        end
        
        function fitted = plotFlowRate(obj)
            fitted = obj.helper_plotFlowRate(2,1);
        end
        
        function fitted = helper_plotFlowRate(obj,lowerlim,ifP)
            fitted = helper_plotFlowRate1(obj,lowerlim,19,1,ifP);
        end
        
        function [fitted, width_b, flowRate_b] = helper_plotFlowRate1(obj,lowerlim,upperlim,fitmode,ifP)
            lower = [0.5,0.4,0];
            upper = [5,1.58,1];
            init = [2,1.1,0.01];
            [fitted, width_b, flowRate_b] = obj.helper_plotFlowRate2(lowerlim,upperlim,fitmode,ifP,upper,lower,init);
        end
        
        function [fitted, width_b, flowRate_b] = helper_plotFlowRate2(obj,lowerlim,upperlim,fitmode,ifP,upper,lower,init)
            flowRate = [];
            error = [];
            width = [];
            for i = 1: length(obj.hopperProperty)
                width = [width, obj.hopperProperty(i).width];
            end
            width = unique(width);
            for w = width
                rate_temp = [];
                std_temp = [];
                for i = 1: length(obj.hopperProperty)
                    if obj.hopperProperty(i).width == w
                        rate_temp = [rate_temp, obj.hopperProperty(i).flow_rate];
                        std_temp = [std_temp, obj.hopperProperty(i).flowStd];
                    end
                end
                flowRate = [flowRate, mean(rate_temp)];
                error = [error, mean(std_temp)];
%                 error = [error, std(rate_temp)/sqrt(length(rate_temp))];

            end
            %12
            filter = width < upperlim & width > lowerlim;
%             filter = width < 11 & width > 4.5;
            width_b = width(filter);
            flowRate_b = flowRate(filter);
            
            if (fitmode <= 2 || fitmode==10)
                filter = flowRate_b > 0.01;
                width_b = width_b(filter);
                flowRate_b = flowRate_b(filter);
            end
            
            if (fitmode==10)
                fitmode=3;
            end
            
            if (fitmode == 1)
            b_fit = 0.5;
             fitfun = fittype( @(a, c, x) abs(c*(x - a).^b_fit));
% % %                 fitfun = fittype( @(a, c, x) c*(x - a).^b_fit);
             fitted = fit( width_b', flowRate_b', fitfun, 'StartPoint', [0.01,1.1]);  
            elseif (fitmode == 2)
                    b_fit = 1.5;
             fitfun = fittype( @(a, c, x) abs(c*(x - a).^b_fit));
% % %                 fitfun = fittype( @(a, c, x) c*(x - a).^b_fit);
             fitted = fit( width_b', flowRate_b', fitfun, 'StartPoint', [0.01,1.1]);  
            elseif (fitmode == 3)
             fo = fitoptions('Method','NonlinearLeastSquares','Upper',upper,'Lower',lower,'StartPoint', init);
             %,'Lower',[0.5,0.5,0],'Upper',[4,1.6,1]
             fitfun = fittype( @(a, b_fit, c, x) abs(c*(x - a).^b_fit),'options',fo);
             fitted = fit( width_b', flowRate_b', fitfun);

            elseif (fitmode == 4)
                    b_fit = 2.5;
             fitfun = fittype( @(a, c, x) abs(c*(x - a).^b_fit));
% % %                 fitfun = fittype( @(a, c, x) c*(x - a).^b_fit);
             fitted = fit( width_b', flowRate_b', fitfun, 'StartPoint', [0.01,1.1]);  
            elseif (fitmode == 5)
              fitfun = fittype( @(a, b_fit, c, x) abs(c*(x - a).^b_fit));
              fitted = fit( width_b', flowRate_b', fitfun, 'StartPoint', [2,0.8,0.008]); 
            end

% % % %              fitfun = fittype( @(a, b_fit, c, x) c*(x - a).^b_fit);
%              fitfun = fittype( @(a, b_fit, c, x) abs(c*(x - a).^b_fit));
%              fitted = fit( width_b', flowRate_b', fitfun, 'StartPoint', [2,0.8,0.008]);
%             fitted = fit( width_b', flowRate_b', fitfun, 'StartPoint', [1,0.5,0.008]);
             disp([fitted])
             if (ifP)
            figure(3);hold on;box on;
            set(gcf,'color','w');
            
            %errorbar(1 - gam/(clog(1,1)*lengthscale(1,1)*coordinate(1,3)),clog_p,error,'o')
            filter2 = width_b > fitted.a;
            width_b = width_b(filter2);
            flowRate_b = flowRate_b(filter2);
            plot(fitted,'k',width_b, flowRate_b);
%             errorbar(width,flowRate,error,'o')
%             plot(fitted,width(width>fitted.a), flowRate(width>fitted.a));
%             scatter(width, flowRate,10,'filled','b');
%             ax = gca;
%             ax.XScale = "log";
%             ax.YScale = "log";
            xlabel("width");
            ylabel("flow rate");
            
            figure(4);hold on;box on;
            set(gcf,'color','w');
            scatter(width_b-fitted.a,abs(flowRate_b/fitted.c));
            ax = gca;
            ax.XScale = "log";
            ax.YScale = "log";
            xlabel("width - a");
            ylabel("flowRate/A");
%             inc = log10(1.4)-0.5*log10(2.5);
%             yy=10^(0.5*log10(9)+inc);
%             plot([2.5,9],[1.4,yy],'Color','k','LineStyle','--')
            inc = log10(3)-1.5*log10(2.5);
            yy=10^(1.5*log10(9)+inc);
            plot([2.5,9],[3,yy],'Color','k','LineStyle','--')
             end
            
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
            eff_width = [];
            for i = 1: length(obj.hopperProperty)
                width = [width, obj.hopperProperty(i).width];
            end
            width = unique(width);
            for w = width
                results = [];
                width_scale = [];
                for i = 1: length(obj.hopperProperty)
                    if obj.hopperProperty(i).width == w
                        results = [results, obj.hopperProperty(i).result];
                        width_scale = [width_scale, obj.hopperProperty(i).eff_width_scale];
                    end
                end
                clog_p = [clog_p, mean(results)];
                eff_width = [eff_width, w * mean(width_scale)];
                error = [error, std(results)/sqrt(size(results, 2))];
            end
            
            %width = eff_width;
            
            fitfun = fittype( @(a, b, x) 1 ./ (1 + exp( (x - a) ./ b)));
            %fitted = fit( (obj.logindex(half:end))', abs(obj.ISF(obj.logindex(half:end))), fitfun, 'StartPoint', [10 * (10 - obj.t_index_j),1]);
            fitted = fit(width', clog_p', fitfun, 'StartPoint', [1,0.1]);
%             fitted.a = 2.6;
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

