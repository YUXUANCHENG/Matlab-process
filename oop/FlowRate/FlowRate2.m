classdef FlowRate2 < FlowRate1
    %FLOWRATE2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        window
    end
    
    methods
        function obj = FlowRate2(count, window)
            obj = obj@FlowRate1(count);
            obj.window = window;           
        end
        
        function rate = help_cal_rate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
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

