classdef Trial_hopper < Trial_DPM
    %TRIAL_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        flow_rate
        gam
        g
        width
        result
        density
        vDistribution
    end
    
    methods
        function obj = Trial_hopper(t_index_i,t_index_j,basefolder,fileReader)
            obj = obj@Trial_DPM(t_index_i,t_index_j,basefolder,fileReader);
            obj.periodic = 0;
        end
        
        function hopperDataProcess(obj)
            obj.fileReader.coordinate(:,1) = obj.fileReader.coordinate(:,1) + obj.lengthscale(end-1);
            obj.lengthscale(end-1) = 2 * obj.lengthscale(end-1);
        end
        
        function readInitial(obj)
            readInitial@Trial_DPM(obj);
            %obj.hopperDataProcess();
        end
        
        function readV0(obj)
            obj.fileReader.readV0();
            obj.kl = obj.fileReader.v0(1);
            obj.gam = obj.fileReader.v0(2);
            obj.g = obj.fileReader.v0(3);
            obj.width = obj.fileReader.v0(4);
            obj.result = obj.fileReader.v0(5);
        end
        
        function readMDdata(obj)
            readMDdata@Trial_DPM(obj);
            %obj.hopperDataProcess();
        end
        
        function createCalculator(obj)
            obj.calculator = Calculator_hopper(obj);
        end
        
        function printCellCount(obj)
            obj.calculator.cell_count();
            temp_c = (obj.calculator.count - obj.calculator.count(1)) * (-1); 
            figure(14), hold on, box on;
            plot(obj.calculator.deltaT, temp_c)
        end
        
        function plotEk(obj)
            Ek = obj.calculator.cal_Ek();
            figure(15), hold on, box on;
            plot(Ek)
            ax = gca;
            %ax.FontSize = 22;
            %ax.XScale = "log";
            ax.YScale = "log";
        end
        
        function flowRate(obj, mode, varargin)
            if isempty(obj.calculator.count)
                obj.calculator.cell_count();
            end
            obj.flow_rate = obj.calculator.flowRate(mode, varargin{1});
        end
        
        function calDensity(obj, timeWindow, heightWindow)
            if isempty(obj.calculator.count)
                obj.calculator.cell_count();
            end
            obj.density = obj.calculator.calDensity(timeWindow, heightWindow);
        end
        
        function calVdistribution(obj,timeWindow, heightWindow)
            if isempty(obj.calculator.count)
                obj.calculator.cell_count();
            end
            obj.vDistribution = obj.calculator.calVdistribution(timeWindow, heightWindow);
        end
        
        function plotDensity(obj)
            figure(13), hold on, box on;
            valueMatrix = obj.density{3};
            valueMatrix = reshape(valueMatrix, [], obj.density{4});
%             h = HeatMap(valueMatrix);
%             %h = heatmap(obj.density{1},obj.density{2},obj.density{3});
%             h.XLabel = '# particle';
%             h.YLabel = 'height';
            colormap('jet')
            imagesc(valueMatrix)
            colorbar
        end
        
        function plotVdistribution(obj)
            figure(7), hold on, box on;
            valueMatrix = obj.vDistribution{3};
            valueMatrix = reshape(valueMatrix, [], obj.vDistribution{4});
%             h = HeatMap(valueMatrix);
%             %h = heatmap(obj.density{1},obj.density{2},obj.density{3});
%             h.XLabel = '# particle';
%             h.YLabel = 'height';
            colormap('jet')
            imagesc(valueMatrix)
            colorbar
        end
        
    end
end

