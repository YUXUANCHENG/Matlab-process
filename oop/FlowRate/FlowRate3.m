classdef FlowRate3 < FlowRate2
    %FLOWRATE3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        originalCount
    end
    
    methods
        
        function rate = cal_rate(obj)
            obj.originalCount = obj.count;
            rate = cal_rate@FlowRate2(obj);
        end
        
        function rate = help_cal_rate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.endtime - obj.starttime <= obj.window
                rate = help_cal_rate@FlowRate1(obj);
            else
                temp_q = [];
                temp_c = [];
                for time_point = obj.starttime : obj.window : obj.endtime - obj.window
                    temp_rate = (obj.count(time_point + obj.window)-obj.count(time_point))/obj.window;
                    if temp_rate > 0
                        temp_q = [temp_q, temp_rate];
                        temp_c = [temp_c, obj.originalCount(time_point)];
                    end
                end
                rate = {temp_q, temp_c};
            end         
        end
        
    end
end

