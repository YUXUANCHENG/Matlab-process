classdef FlowRate1 < handle
    %FLOWRATE1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        count
        starttime
        endtime
    end
    
    methods
        
        function obj = FlowRate1(count)
            %FLOWRATE2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.count = count;
        end
        
        function rate = cal_rate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.count =  (obj.count - obj.count(1)) * (-1); 
            
            obj.starttime = find(obj.count > 0);
            if size(obj.starttime,1) > 0
                obj.starttime = obj.starttime(1);
                obj.endtime = find(obj.count == obj.count(end));
                obj.endtime = obj.endtime(1);
                if obj.endtime > obj.starttime
                    rate = help_cal_rate(obj);
                else
                    rate = 0;
                end
            else
                rate = 0;
            end
        end
        
        function rate = help_cal_rate(obj)
            rate = (obj.count(end)-obj.count(obj.starttime))/(obj.endtime - obj.starttime);
        end
    end
end

