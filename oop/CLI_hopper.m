classdef CLI_hopper < CLI
    %CLI_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hopperProperty = []
    end
    
    methods
        function calHopperProperty(obj,index_i, index_j, mode, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for t_index_i =0: index_i
                for t_index_j = 0: index_j
                    try
                        fileReader = FileReader();
                        trial = Trial_hopper(t_index_i,t_index_j,obj.basefolder,fileReader);
                        %trial.plotInitial();
                        trial.readV0();
                        trial.readMDdata();
                        trial.plotLastFrame(2);
                        trial.createCalculator();
                        trial.printCellCount();
                        trial.flowRate(mode,varargin);
                        %obj.clog = [obj.clog;[trial.fileReader.v0, trial.flow_rate]];
                        obj.hopperProperty = [obj.hopperProperty, trial];
                    catch e
                        fprintf(1,"%s", e.message);
                        delete(fileReader);
                        continue
                    end
                    delete(fileReader);
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
        
    end
end

