classdef Calculator_hopper < Calculator
    %CALCULATOR_HOPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        count
        aboveOpenning
        deltaT
        flowRateCalculator
        density
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
        
        function angle = cal_vel_distri_in_space(obj,direction)
            tlim = obj.trial.frames;
            figure(1);hold on;box on;
            set(gcf,'color','w');
            space = 0.5;
%             space = 2;
            if direction == 1
%                 maxIter = 9;
%                 maxIter = 4;
                  maxIter = 19;
            else
                maxIter = 0;
            end
            gridV = [];
            
            for i = 0:maxIter
                index = 0;
                test = 0;
%                 for t = 4:19
                for t = round(tlim/3):tlim-0
%                   for t = tlim-5:tlim-0
                    index = index + 1;
                    [returnX, returnY] = obj.help_cal_vel_distri_in_space2(direction,t,i,space);
                    test = test + returnY;
                    
                end
                test = test/index;
                gridV = [gridV test];
%                 plot(returnX,test)
            end
            figure(3);
            maxV = max(gridV,[],'all');
            dataV=gridV;
%             dataV = (abs(gridV./maxV));
%             dataV = log10(abs(gridV./maxV));
%             threash = -0.9;
%             threash = -1.1;
%               threash = 0.05;
%              threash = 0.1;
%             dataV(dataV<threash) = threash;
%             heatmap(dataV,'Colormap', jet);
%             heatmap(dataV);
            contour(dataV);
            grid off;
%             heatmap(log10(abs(gridV)),'Colormap', jet);
            figure(5);set(gcf,'color','w');
            [px,py] = gradient(dataV);
            offset = 1;
            px=px(:,offset:end);py=py(:,offset:end);dataV=dataV(:,offset:end);
            contour(dataV,0.1:0.1:1)
            hold on
            quiver(px,py)
            hold off
            
%             %%angle based on gradiant
%             magP = sqrt(px.^2+py.^2);
%             tanAngle = rad2deg(atan(py./px));
%             width = obj.trial.width;
%             widthInGrid = ceil((1*width*(1+1.4)/4)/space)+1;
%             upperlim = 15 + widthInGrid;
%             lowerlim = 15 - widthInGrid -1;
%             upperRegionSum = sum(magP(upperlim:end,:),'all');
%             lowerRegionSum = sum(magP(1:lowerlim,:),'all');
% %             avgTanUpper = sum(tanAngle(upperlim:end,:).* magP(upperlim:end,:)/upperRegionSum,'all');
% %             avgTanLower = sum(tanAngle(1:lowerlim,:).* magP(1:lowerlim,:)/lowerRegionSum,'all');
%             avgTanUpper = rad2deg(atan(sum(py(upperlim:end,:),'all')/sum(px(upperlim:end,:),'all')));
%             avgTanLower = rad2deg(atan(sum(py(1:lowerlim,:),'all')/sum(px(1:lowerlim,:),'all')));
%             uperAngle = (-avgTanUpper);
%             lowerAngle = (avgTanLower);
%             if (uperAngle<0) 
%                 uperAngle = lowerAngle;
%             end
%             if (lowerAngle<0) 
%                 lowerAngle = uperAngle;
%             end
%             angle = 90-(uperAngle+lowerAngle)/2;
            
%             %%angle based on threashold
%             boolDataV = dataV>threash;
%             upper=[];
%             lower=[];
%             for xIndex = 1: size(dataV,2)
%                 tt = find(boolDataV(:,xIndex));
%                 maxValue = max(tt);
%                 minValue = min(tt);
%                 if (maxValue == minValue)
%                     break
%                 end
%                 if (length(tt)<2)
%                     break
%                 end
%                 upper = [upper, max(tt)];
%                 lower = [lower, min(tt)];
%             end
%             fit1 = polyfit(1: length(upper),upper,1);
%             fit2 = polyfit(1: length(lower),lower,1);
%             if (fit2(1)>0)
%                 fit2 = -fit1;
%             end
%             meanSlop = (fit1(1)-fit2(1))/2;
%             %disp(atan(meanSlop)/3.1415 * 180)
%             angle = atan(meanSlop)/3.1415 * 180;
        end
        
        function [returnX, returnY] = help_cal_vel_distri_in_space1(obj,direction,time,i,space)
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end
            cEnd = time * obj.trial.Ncell;
            cStart = cEnd - obj.trial.Ncell + 1;
            calA = obj.trial.fileReader.cal_A(cStart:cEnd);
            velU = obj.trial.fileReader.vel(cStart:cEnd, 1);
            velV = obj.trial.fileReader.vel(cStart:cEnd, 2);
            
%             filter = obj.x_comp(time,:) < obj.trial.lengthscale(end-1);
%             meanU = mean(velU(filter),'all');
%             meanV = mean(velV(filter),'all');
            velocity = sqrt(velU.^2 + velV.^2);
%             figure(1);hold on;box on;
%             set(gcf,'color','w');
            
            if direction == 0
            space_tick = 0:-3:-21;
            vel = zeros(length(space_tick)-1,1);
            calAseg = zeros(length(space_tick)-1,1);
            num_in_bin = zeros(length(space_tick)-1,1);
            real_width = obj.trial.width * 2 / (1 + 1.4);
            upper = (obj.trial.lengthscale(end) + real_width) / 2;
            lower = (obj.trial.lengthscale(end) - real_width) / 2;
            for index = 1 : length(obj.x_comp(end,:))
                if obj.x_comp(time,index) < obj.trial.lengthscale(end-1) && obj.y_comp(time,index) < upper && obj.y_comp(time,index) > lower
                    for bin_index = 1: length(space_tick)-1
                        if obj.x_comp(time,index) < space_tick(bin_index) && obj.x_comp(time,index) > space_tick(bin_index+1)
                            num_in_bin(bin_index) = num_in_bin(bin_index) + 1;
%                             vel(bin_index) = vel(bin_index) + velU(index);
%                             vel(bin_index) = vel(bin_index) + velV(index);
%                             vel(bin_index) = vel(bin_index) + sqrt(velV(index).^2+velU(index).^2);
%                             vel(bin_index) = vel(bin_index) + (velV(index)-meanV).^2 + (velU(index)-meanU).^2;
                            calAseg(bin_index) = calAseg(bin_index) + calA(index);
                        end
                    end
                end
            end
            vel = vel./num_in_bin;
            calAseg = calAseg./num_in_bin;
            returnX = space_tick(1:end-1);
            returnY = vel;
%             plot(space_tick(1:end-1),vel)
            else
%             for i = 0:0
%             for i = 0:0
            
            space_tick = 0:space:obj.trial.lengthscale(end);
            vel = zeros(length(space_tick)-1,1);
            meanU = zeros(length(space_tick)-1,1);
            meanV = zeros(length(space_tick)-1,1);
            calAseg = zeros(length(space_tick)-1,1);
            num_in_bin = zeros(length(space_tick)-1,1);
            upper =  - i*space;
            lower = -space*(i+1);
            for index = 1 : length(obj.x_comp(end,:))
                if obj.x_comp(time,index) < obj.trial.lengthscale(end-1) && obj.x_comp(time,index) < upper && obj.x_comp(time,index) > lower
                    for bin_index = 1: length(space_tick)-1
                        if obj.y_comp(time,index) > space_tick(bin_index) && obj.y_comp(time,index) < space_tick(bin_index+1)
                            num_in_bin(bin_index) = num_in_bin(bin_index) + 1;
                            meanU(bin_index) = meanU(bin_index) + velU(index);
                            meanV(bin_index) = meanV(bin_index) + velV(index);
%                             vel(bin_index) = vel(bin_index) + velU(index);
%                             vel(bin_index) = vel(bin_index) + velV(index);
%                             vel(bin_index) = vel(bin_index) + (velV(index)-meanV).^2 + (velU(index)-meanU).^2;
%                             calAseg(bin_index) = calAseg(bin_index) + calA(index);
                        end
                    end
                end
            end
            meanU = meanU./num_in_bin;meanU(isnan(meanU)) = 0;
            meanV = meanV./num_in_bin;meanV(isnan(meanV)) = 0;
            num_in_bin = zeros(length(space_tick)-1,1);
            
            for index = 1 : length(obj.x_comp(end,:))
                if obj.x_comp(time,index) < obj.trial.lengthscale(end-1) && obj.x_comp(time,index) < upper && obj.x_comp(time,index) > lower
                    for bin_index = 1: length(space_tick)-1
                        if obj.y_comp(time,index) > space_tick(bin_index) && obj.y_comp(time,index) < space_tick(bin_index+1)
                            num_in_bin(bin_index) = num_in_bin(bin_index) + 1;
%                             vel(bin_index) = vel(bin_index) + velU(index);
%                             vel(bin_index) = vel(bin_index) + velV(index);
                            vel(bin_index) = vel(bin_index) + (velV(index)-meanV(bin_index)).^2 + (velU(index)-meanU(bin_index)).^2;
                            calAseg(bin_index) = calAseg(bin_index) + calA(index);
                        end
                    end
                end
            end
  
            vel = vel./num_in_bin;
            calAseg = calAseg./num_in_bin;
            vel(isnan(vel)) = 0;
       
            returnX = space_tick(1:end-1)-obj.trial.lengthscale(end)/2;
            returnY = vel;
%             plot(space_tick(1:end-1)-obj.trial.lengthscale(end)/2,calAseg)
%             plot(space_tick(1:end-1)-obj.trial.lengthscale(end)/2,vel)
%             end
            
            end
        end
       function [returnX, returnY] = help_cal_vel_distri_in_space2(obj,direction,time,i,space)
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end
            cEnd = time * obj.trial.Ncell;
            cStart = cEnd - obj.trial.Ncell + 1;
            calA = obj.trial.fileReader.cal_A(cStart:cEnd);
            velU = obj.trial.fileReader.vel(cStart:cEnd, 1);
            velV = obj.trial.fileReader.vel(cStart:cEnd, 2);
            
%             filter = obj.x_comp(time,:) < obj.trial.lengthscale(end-1);
%             meanU = mean(velU(filter),'all');
%             meanV = mean(velV(filter),'all');
            velocity = sqrt(velU.^2 + velV.^2);
%             figure(1);hold on;box on;
%             set(gcf,'color','w');
            
            if direction == 0
                jj=1;
%             space_tick = 0:-3:-21;
%             vel = zeros(length(space_tick)-1,1);
%             calAseg = zeros(length(space_tick)-1,1);
%             num_in_bin = zeros(length(space_tick)-1,1);
%             real_width = obj.trial.width * 2 / (1 + 1.4);
%             upper = (obj.trial.lengthscale(end) + real_width) / 2;
%             lower = (obj.trial.lengthscale(end) - real_width) / 2;
%             for index = 1 : length(obj.x_comp(end,:))
%                 if obj.x_comp(time,index) < obj.trial.lengthscale(end-1) && obj.y_comp(time,index) < upper && obj.y_comp(time,index) > lower
%                     for bin_index = 1: length(space_tick)-1
%                         if obj.x_comp(time,index) < space_tick(bin_index) && obj.x_comp(time,index) > space_tick(bin_index+1)
%                             num_in_bin(bin_index) = num_in_bin(bin_index) + 1;
% %                             vel(bin_index) = vel(bin_index) + velU(index);
% %                             vel(bin_index) = vel(bin_index) + velV(index);
% %                             vel(bin_index) = vel(bin_index) + sqrt(velV(index).^2+velU(index).^2);
% %                             vel(bin_index) = vel(bin_index) + (velV(index)-meanV).^2 + (velU(index)-meanU).^2;
%                             calAseg(bin_index) = calAseg(bin_index) + calA(index);
%                         end
%                     end
%                 end
%             end
%             vel = vel./num_in_bin;
%             calAseg = calAseg./num_in_bin;
%             returnX = space_tick(1:end-1);
%             returnY = vel;
%             plot(space_tick(1:end-1),vel)
            else
            space_tick = 0:space:obj.trial.lengthscale(end);
            rho = zeros(length(space_tick)-1,1);
            px = zeros(length(space_tick)-1,1);
            py = zeros(length(space_tick)-1,1);
            VVx = zeros(length(space_tick)-1,1);
            VVy = zeros(length(space_tick)-1,1);
            kStress = zeros(length(space_tick)-1,1);
            vel = zeros(length(space_tick)-1,1);
            meanU = zeros(length(space_tick)-1,1);
            meanV = zeros(length(space_tick)-1,1);
            calAseg = zeros(length(space_tick)-1,1);
            num_in_bin = zeros(length(space_tick)-1,1);
            limUp = obj.trial.lengthscale(end-1)-0.24;
            upper =  limUp- i*space;
            lower = limUp -space*(i+1);
%             
            for bin_index = 1: length(space_tick)-1
                yloc = space_tick(bin_index);
                xloc = upper;
                for index = 1 : length(obj.x_comp(end,:))
                    rr = sqrt((obj.y_comp(time,index)-yloc).^2 + (obj.x_comp(time,index)-xloc).^2);
                    cgCoeff = (1/(sqrt(2*pi)*0.6)^3)*exp(-rr^2/(2*0.6^2));
                    rho(bin_index) = rho(bin_index) + cgCoeff;
                    px(bin_index) = px(bin_index) + velU(index)*cgCoeff;
                    py(bin_index) = py(bin_index) + velV(index)*cgCoeff;
                end
            end

            VVx = px./rho; VVx(isnan(VVx)) = 0;
            VVy = py./rho; VVy(isnan(VVy)) = 0;
            
            for bin_index = 1: length(space_tick)-1
                yloc = space_tick(bin_index);
                xloc = upper;
                for index = 1 : length(obj.x_comp(end,:))
                    rr = sqrt((obj.y_comp(time,index)-yloc).^2 + (obj.x_comp(time,index)-xloc).^2);
                    cgCoeff = (1/(sqrt(2*pi)*0.6)^3)*exp(-rr^2/(2*0.6^2));
%                     kStress(bin_index) = kStress(bin_index) + cgCoeff *(((velV(index)-VVy(bin_index))/VVy(bin_index)).^2 + ((velU(index)-VVx(bin_index))/VVx(bin_index)).^2);
                    kStress(bin_index) = kStress(bin_index) + cgCoeff *(((velV(index)-VVy(bin_index))).^2 + ((velU(index)-VVx(bin_index))).^2);
                
                end
            end
  
%             vel = vel./num_in_bin;
%             calAseg = calAseg./num_in_bin;
%             vel(isnan(vel)) = 0;
       
            returnX = space_tick(1:end-1)-obj.trial.lengthscale(end)/2;
%             returnY = px;
            returnY = kStress;



            
            end
        end

          function [returnX, returnY] = help_cal_vel_distri_in_space(obj,direction,time,i,space)
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end
            cEnd = time * obj.trial.Ncell;
            cStart = cEnd - obj.trial.Ncell + 1;
            calA = obj.trial.fileReader.cal_A(cStart:cEnd);
            velU = obj.trial.fileReader.vel(cStart:cEnd, 1);
            velV = obj.trial.fileReader.vel(cStart:cEnd, 2);
            velocity = sqrt(velU.^2 + velV.^2);
            
            start_point = 1 + obj.trial.N * ( time - 1 );
            end_point = obj.trial.N * time;
            x_f = obj.trial.fileReader.coordinate(start_point:end_point,1);
            y_f = obj.trial.fileReader.coordinate(start_point:end_point,2);
%             figure(1);hold on;box on;
%             set(gcf,'color','w');
            
            space=0.5;
            space_tick = 0:space:obj.trial.lengthscale(end);
            vel = zeros(length(space_tick)-1,1);
            meanU = zeros(length(space_tick)-1,1);
            meanV = zeros(length(space_tick)-1,1);
            calAseg = zeros(length(space_tick)-1,1);
            num_in_bin = zeros(length(space_tick)-1,1);
            upper = 0 - i*space;
            lower = -space*(i+1);
            for ci = 1 : length(obj.x_comp(end,:))
                index = sum(obj.trial.fileReader.lengthscale(1:ci),'all');
                start_point_last = 1 + index - obj.trial.fileReader.lengthscale(ci);
                end_point_last = index;
                x = x_f(start_point_last:end_point_last);
                y = y_f(start_point_last:end_point_last);
                x = [x; x(1)];
                y = [y; y(1)];
                
                     boxY = (space_tick(1:end-1)+space_tick(2:end))/2;
                     boxX = ones(length(space_tick)-1,1)*(upper+lower)/2;
                     ifinside = inpolygon(boxX,boxY,x,y);
               
                        num_in_bin(ifinside) = num_in_bin(ifinside) + 1;
                        meanU(ifinside) = meanU(ifinside) + velU(ci);
                        meanV(ifinside) = meanV(ifinside) + velV(ci);
%                         vel(ifinside) = vel(ifinside) + velU(ci);
    %                     vel(bin_index) = vel(bin_index) + velV(index);
                        %calAseg(ifinside) = calAseg(ifinside) + calA(ci);
                
           
             
            end
            meanU = meanU./num_in_bin;meanU(isnan(meanU)) = 0;
            meanV = meanV./num_in_bin;meanV(isnan(meanV)) = 0;
            num_in_bin = zeros(length(space_tick)-1,1);           
            for ci = 1 : length(obj.x_comp(end,:))
                index = sum(obj.trial.fileReader.lengthscale(1:ci),'all');
                start_point_last = 1 + index - obj.trial.fileReader.lengthscale(ci);
                end_point_last = index;
                x = x_f(start_point_last:end_point_last);
                y = y_f(start_point_last:end_point_last);
                x = [x; x(1)];
                y = [y; y(1)];
                
                     boxY = (space_tick(1:end-1)+space_tick(2:end))/2;
                     boxX = ones(length(space_tick)-1,1)*(upper+lower)/2;
                     ifinside = inpolygon(boxX,boxY,x,y);
               
                        num_in_bin(ifinside) = num_in_bin(ifinside) + 1;
                        vel(ifinside) = vel(ifinside) + (velV(ci)-meanV(ifinside)).^2 + (velU(ci)-meanU(ifinside)).^2;
    %                     vel(bin_index) = vel(bin_index) + velV(index);
                        %calAseg(ifinside) = calAseg(ifinside) + calA(ci);
                
           
             
            end
  
            vel = vel./num_in_bin;
            calAseg = calAseg./num_in_bin;
            vel(~(vel>0)) = 1e-5;
       
            returnX = space_tick(1:end-1)-obj.trial.lengthscale(end)/2;
            returnY = vel;
%             plot(space_tick(1:end-1)-obj.trial.lengthscale(end)/2,calAseg)
%             plot(space_tick(1:end-1)-obj.trial.lengthscale(end)/2,vel)
%             end
            
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
        
        function density = calDensity(obj, timeWindow, heightWindow)
            density = {};
            numberInHopper = [];
            height = [];
            densittAtHeight = [];
            count1 = 0;
            for ii = 1: timeWindow: length(obj.count)
                %if (size(numberInHopper,2) > 0) && (obj.count(ii) == numberInHopper(end))
                if ismember(obj.count(ii), numberInHopper) || obj.count(ii) == obj.count(1)
                    continue
                end
                
                for xx = 0 : -heightWindow : -obj.trial.lengthscale(end) * 3
                    % x axis, # particles in hopper
                    numberInHopper= [numberInHopper, obj.count(ii)];
                    height = [height, xx];
                    num = sum(obj.x_comp(ii,:) < xx & obj.x_comp(ii,:) > xx - heightWindow);
                    densittAtHeight= [densittAtHeight, num/(obj.trial.lengthscale(end)*heightWindow)];
                end
                count1 = count1 + 1;
            end
            density{1} = numberInHopper;
            density{2} = height;
            density{3} = densittAtHeight;  
            density{4} = count1;
        end
        
        function density = calDensityAtOpening(obj, timeWindow, heightWindow)
            density = {};
            numberInHopper = [];
            height = [];
            densittAtHeight = [];
            count1 = 0;
            for ii = 3  
                for xx = 0
                    % x axis, # particles in hopper
                    numberInHopper= [numberInHopper, obj.count(ii)];
                    height = [height, xx];
                    num = sum(obj.x_comp(ii,:) < xx & obj.x_comp(ii,:) > xx - heightWindow);
                    densittAtHeight= [densittAtHeight, num/(obj.trial.lengthscale(end)*heightWindow)];
                end
                count1 = count1 + 1;
            end
            density{1} = numberInHopper;
            density{2} = height;
            density{3} = densittAtHeight;  
            density{4} = count1;
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

