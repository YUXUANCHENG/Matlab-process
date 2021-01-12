classdef Calculator_hopper < Calculator
    %CALCULATOR_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        count
        deltaT
        flowRateCalculator
    end
    
    methods
        function cell_count(obj)
            xcomp=zeros(obj.trial.frames,obj.trial.Ncell);
            ycomp=zeros(obj.trial.frames,obj.trial.Ncell);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.N * ( i - 1 );
                end_point = obj.trial.N * i;
                [xcomp_t,ycomp_t] = obj.cal_c_pos(obj.trial.fileReader.coordinate(start_point:end_point,1),obj.trial.fileReader.coordinate(start_point:end_point,2));
                xcomp(i,:)= xcomp_t;
                ycomp(i,:)= ycomp_t;
            end
            % create MSD array (y-axis of MSD plot)
            obj.count = zeros(obj.trial.frames,1);
            NT = length(obj.count);
            % logindex = unique(round(logspace(0, log10(NT),200)));
            % loop over the different possible time windows, calculate MSD for each
            % time window size
            for ii = 1:NT
                % store in MSD array
                obj.count(ii) = sum(xcomp(ii,:) < obj.trial.lengthscale(end-1));
            end

            % create deltaT array, using a for loop or vectorization
            obj.deltaT = 1:NT;
            % obj.count =  (obj.count - obj.count(1)) * (-1); 
        end
        
        function rate = flowRate(obj, mode, varargin)
            if mode == 1
                obj.flowRateCalculator = FlowRate1(obj.count);
            elseif mode == 2
                obj.flowRateCalculator = FlowRate2(obj.count, varargin{1}{1});
            elseif mode == 3
                obj.flowRateCalculator = FlowRate3(obj.count, varargin{1}{1});
            end
            rate = obj.flowRateCalculator.cal_rate();
        end
        
        function Ek = cal_Ek(obj)
            Ek = zeros(obj.trial.frames,1);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.Ncell * ( i - 1 );
                end_point = obj.trial.Ncell * i;
                Ek(i) =  sum(obj.trial.fileReader.vel(start_point:end_point),'all');
            end   
        end
        
    end
end

