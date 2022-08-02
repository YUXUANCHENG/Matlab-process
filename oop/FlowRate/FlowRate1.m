classdef FlowRate1 < handle
    %FLOWRATE1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        calculator
        count
        starttime
        endtime
        cutoff = 0.9
%         cutoff = 1
    end
    
    methods
        
        function obj = FlowRate1(calculator)
            %FLOWRATE2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.calculator = calculator;
            obj.count = obj.calculator.count;
        end
        
        function rate = cal_rate(obj)
            % function for calculation flow rate
            obj.count =  (obj.count - obj.count(1)) * (-1); 
            
            obj.starttime = find(obj.count > 0);
            if size(obj.starttime,1) > 0
                obj.starttime = obj.starttime(1);
                obj.endtime = find(obj.count == floor(obj.count(end)*obj.cutoff));
                try
                    obj.endtime = obj.endtime(1);
                    %obj.endtime = obj.endtime(end);
                catch e
                    obj.endtime = 0;
                end
                if obj.endtime > obj.starttime
                    % call the help function, this will be overide in
                    % different implementations
                    rate = help_cal_rate(obj);
                else
                    rate = 0;
                end
            else
                rate = 0;
            end
        end
        
        function rate = help_cal_rate(obj)
            % specific version of cal_rate
            rate = (floor(obj.count(end)*obj.cutoff)-obj.count(obj.starttime))/(obj.endtime - obj.starttime);
        end
    end
end

