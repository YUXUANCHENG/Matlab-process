classdef Calculator_hopper < Calculator
    %CALCULATOR_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        count
        aboveOpenning
        deltaT
        flowRateCalculator
    end
    
    methods
        function cell_count(obj)
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end
            % create MSD array (y-axis of MSD plot)
            obj.count = zeros(obj.trial.frames,1);
            NT = length(obj.count);
            % logindex = unique(round(logspace(0, log10(NT),200)));
            % loop over the different possible time windows, calculate MSD for each
            % time window size
            for ii = 1:NT
                % store in MSD array
                obj.count(ii) = sum(obj.x_comp(ii,:) < obj.trial.lengthscale(end-1));
            end
            % create deltaT array, using a for loop or vectorization
            obj.deltaT = 1:NT;
            % obj.count =  (obj.count - obj.count(1)) * (-1); 
        end
        
        function cal_above_openning(obj)
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end
            % create MSD array (y-axis of MSD plot)
            obj.aboveOpenning = zeros(obj.trial.frames,1);
            NT = length(obj.count);
            
            real_width = obj.trial.width * 2 / (1 + 1.4);
            upper = (obj.trial.lengthscale(end) + real_width) / 2;
            lower = (obj.trial.lengthscale(end) - real_width) / 2;
            
            for ii = 1: NT
                for index = 1 : length(obj.x_comp(ii,:))
                    if obj.x_comp(ii,index) < obj.trial.lengthscale(end-1) && obj.y_comp(ii,index) < upper && obj.y_comp(ii,index) > lower
                        obj.aboveOpenning(ii) = obj.aboveOpenning(ii) + 1;
                    end
                end
            end       
        end
        
        function rate = flowRate(obj, mode, varargin)
            if mode == 1
                obj.flowRateCalculator = FlowRate1(obj);
            elseif mode == 2
                obj.flowRateCalculator = FlowRate2(obj, varargin{1}{1});
            elseif mode == 3
                obj.flowRateCalculator = FlowRate3(obj, varargin{1}{1});
            elseif mode == 4
                obj.cal_above_openning();
                obj.flowRateCalculator = FlowRate4(obj, varargin{1}{1});
            end
            rate = obj.flowRateCalculator.cal_rate();
        end
        
        function Ek = cal_Ek(obj)
            Ek = zeros(obj.trial.frames,1);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.Ncell * ( i - 1 );
                end_point = obj.trial.Ncell * i;
                vx = obj.trial.fileReader.vel(start_point:end_point, 1);
                vy = obj.trial.fileReader.vel(start_point:end_point, 2);
                Ek(i) =  sum(vx.^2 + vy.^2,'all');
            end   
        end
        
    end
end

