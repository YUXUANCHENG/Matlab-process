classdef Trial_hopper < Trial_DPM
    %TRIAL_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        flow_rate
        gam
        g
        width
        eff_width_scale = 1
        result
        density
        flow
        flowStd
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
            obj.width = roundn(obj.fileReader.v0(4),-4);
            obj.result = obj.fileReader.v0(5);
        end
        
        function readAndPlotFriction(obj)
            friction = obj.fileReader.readFriction();
            figure(5);hold on;box on;
            set(gcf,'color','w');
            plot(friction(:,1),friction(:,2),'--','LineWidth',3)
            xlabel("force");
            ylabel("Ek");
        end
        
        function readFlowR(obj)
            [obj.flow, obj.flow_rate] = obj.fileReader.readFlowR();
            obj.flowStd = std(obj.flow);
        end
        
        function plotFlow(obj)
            figure(1); box on;
%             gcf.Position(3:4)=[550,100];
            set(gcf,'color','w');
            set(gcf,'position',[100,100,500,100])
            hasFlow = obj.flow > 0;
            hasFlow = int8(hasFlow)';
            h = heatmap(hasFlow)
            grid(h, 'off')
%             histogram(hasFlow)
            Ax = gca;
            Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
            Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
        end
        
        function readMDdata(obj)
            readMDdata@Trial_DPM(obj);
            %obj.hopperDataProcess();
        end
        
        function effR = calEffR(obj)
            if (isempty(obj.fileReader.coordinate))
                obj.readMDdata();
            end
            obj.createCalculator();
            if obj.result
                i = floor(obj.frames);
            else
                obj.calculator.cell_count();
                [v, i] = min(abs(obj.calculator.count - obj.Ncell* 2/3));
            end
            start_point = 1 + obj.N * ( i - 1 );
            end_point = obj.N * i;
            x_pos = obj.fileReader.coordinate(start_point:end_point,1);
            y_pos = obj.fileReader.coordinate(start_point:end_point,2);
            [xcomp,ycomp] = obj.calculator.help_cal_c_pos(x_pos, y_pos);
            areas = obj.calculator.cal_area(x_pos, y_pos);
%             smallestR = obj.calculator.calSmallestR(x_pos, y_pos);
            longAxis = obj.calculator.calLongAxis(x_pos, y_pos);
            
            index = obj.lengthscale(end-1)/2  < xcomp & xcomp < obj.lengthscale(end-1);
            longAxis = longAxis(index);
            areas = areas(index);
            xcomp = xcomp(index);
            [temp, selected] = maxk(xcomp,4);
            %area_filter = areas > 0;
            %avg_area = mean(areas(index & area_filter),'all');
            %effR = sqrt(avg_area/pi);
            %areas = areas(index & area_filter);
            
%             Rfilter = smallestR(index);
            %effR = mean(sqrt(areas./pi),'all');
%             effR = mean(Rfilter,'all');
            effR = mean(areas(selected)/longAxis(selected),'all')/2;
            
        end
        
        function cal_eff_width_scale(obj)
            effR = obj.calEffR();
            obj.eff_width_scale = ((1 + 1.4)/(2 * 2))/effR;
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
            set(gcf,'color','w');
            plot(Ek)
            ax = gca;
            %ax.FontSize = 22;
            %ax.XScale = "log";
            ax.YScale = "log";
            xlabel("time");
            ylabel("Ek");
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
        
        function volumnRate(obj)
            if isempty(obj.calculator.count)
                obj.calculator.cell_count();
            end
            heightWindow = 3;
            obj.density = obj.calculator.calDensityAtOpening(1, heightWindow);
            scale = 1/(obj.density{3}(1)* (pi*(0.5^2+0.7^2)/2));
            obj.flow_rate = obj.flow_rate * scale;
        end
        
    end
end

