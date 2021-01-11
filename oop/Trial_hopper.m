classdef Trial_hopper < Trial
    %TRIAL_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        flow_rate
    end
    
    methods
        function obj = Trial_hopper(t_index_i,t_index_j,basefolder,fileReader)
            obj = obj@Trial(t_index_i,t_index_j,basefolder,fileReader);
            obj.periodic = 0;
        end
        
        function hopperDataProcess(obj)
            obj.fileReader.coordinate(:,1) = obj.fileReader.coordinate(:,1) + obj.lengthscale(end-1);
            obj.lengthscale(end-1) = 2 * obj.lengthscale(end-1);
        end
        
        function readInitial(obj)
            readInitial@Trial(obj);
            obj.hopperDataProcess();
        end
        
        function readMDdata(obj)
            readMDdata@Trial(obj);
            obj.hopperDataProcess();
        end
        
        function createCalculator(obj)
            obj.calculator = Calculator_hopper(obj);
        end
        
        function printCellCount(obj)
            obj.calculator.cell_count();
            figure(14), hold on, box on;
            plot(obj.calculator.deltaT, obj.calculator.count)
        end
        
        function flowRate(obj, mode, varargin)
            if isempty(obj.calculator.count)
                obj.calculator.cell_count();
            end
            obj.flow_rate = obj.calculator.flowRate(mode, varargin{1});
        end
        
    end
end

