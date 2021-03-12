classdef FlowRate2 < FlowRate1
    %flow rate with specific window size
    
    properties
        window
    end
    
    methods
        function obj = FlowRate2(calculator, window)
            obj = obj@FlowRate1(calculator);
            obj.window = window;           
        end
        
        function rate = help_cal_rate(obj)
            %overide help function
            if obj.endtime - obj.starttime <= obj.window
                rate = help_cal_rate@FlowRate1(obj);
            else
                temp_q = [];
                for time_point = obj.starttime : obj.window : obj.endtime - obj.window
                    temp_rate = (obj.count(time_point + obj.window)-obj.count(time_point))/obj.window;
                    if temp_rate > 0
                        temp_q = [temp_q, temp_rate];
                    end
                end
                rate = mean(temp_q,'all');
            end         
        end
    end
end

