classdef FlowRate4 < FlowRate3
    %flow rate with particles above opening
    
    properties

    end
    
    methods
        function rate = help_cal_rate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.endtime - obj.starttime <= obj.window
                rate = {help_cal_rate@FlowRate1(obj), obj.calculator.aboveOpenning(end)};
            else
                temp_q = [];
                temp_c = [];
                for time_point = obj.starttime : obj.window : obj.endtime - obj.window
                    temp_rate = (obj.count(time_point + obj.window)-obj.count(time_point))/obj.window;
                    if temp_rate >= 0
                        temp_q = [temp_q, temp_rate];
                        temp_c = [temp_c, obj.calculator.aboveOpenning(time_point + floor(obj.window/2))];
                    end
                end
                rate = {temp_q, temp_c};
            end         
        end
    end
end

