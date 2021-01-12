classdef Calculator < handle
    %CALCULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        trial
    end
    
    methods
        function obj = Calculator(trial)
            %CALCULATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.trial = trial;
        end
        
        function [MSD,deltaT] = cal_msd(obj)

            xcomp = zeros(obj.trial.frames, obj.trial.Ncell);
            ycomp = zeros(obj.trial.frames, obj.trial.Ncell);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.N * ( i - 1 );
                end_point = obj.trial.N * i;
                [xcomp_t,ycomp_t] = obj.cal_c_pos(obj.trial.fileReader.coordinate(start_point:end_point,1),obj.trial.fileReader.coordinate(start_point:end_point,2));
                xcomp(i,:)= xcomp_t;
                ycomp(i,:)= ycomp_t;
            end

            % create MSD array (y-axis of MSD plot)
            MSD = zeros(round(9 * obj.trial.frames/10),1);
            NT = length(MSD);
            logindex = unique(round(logspace(0, log10(NT),1000)));
            % loop over the different possible time windows, calculate MSD for each
            % time window size
            for ii = logindex
                % calculate x displacements, separated by ii indices
                dx = xcomp(1+ii:end,:) - xcomp(1:end-ii,:);
                % calculate y displacements similarly
                dy = ycomp(1+ii:end,:) - ycomp(1:end-ii,:);
                % take mean over all displacements
                dispMean = mean(dx.^2 + dy.^2,'all');
                % store in MSD array
                MSD(ii) = dispMean;
            end

            % create deltaT array, using a for loop or vectorization
            deltaT = 1:NT;
        end
        
        function [xcomp,ycomp] = cal_c_pos(obj, xpos_at_frame, ypos_at_frame)
            xcomp = zeros(1, obj.trial.Ncell);
            ycomp = zeros(1, obj.trial.Ncell);
            for ci = 1: obj.trial.Ncell
                index = sum(obj.trial.lengthscale(1:ci),'all');
                start_point_last = 1 + index - obj.trial.lengthscale(ci);
                end_point_last = index;
                cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
                cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
                xcomp(ci) = cx_tmp;
                ycomp(ci) = cy_tmp;
            end
        end
        
    end
end

