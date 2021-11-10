classdef Calculator < handle
    %CALCULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        trial
        x_comp
        y_comp
        x_v
    end
    
    methods
        function obj = Calculator(trial)
            %CALCULATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.trial = trial;
        end
        
        function [MSD,deltaT] = cal_msd(obj)
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end

            % create MSD array (y-axis of MSD plot)
            MSD = zeros(round(9 * obj.trial.frames/10),1);
            NT = length(MSD);
            logindex = unique(round(logspace(0, log10(NT),1000)));
            % loop over the different possible time windows, calculate MSD for each
            % time window size
            for ii = logindex
                % calculate x displacements, separated by ii indices
                dx = obj.x_comp(obj.trial.offset+ii:end,:) - obj.x_comp(obj.trial.offset:end-ii,:);
                % calculate y displacements similarly
                dy = obj.y_comp(obj.trial.offset+ii:end,:) - obj.y_comp(obj.trial.offset:end-ii,:);
                % take mean over all displacements
                dispMean = mean(dx.^2 + dy.^2,'all');
                % store in MSD array
                MSD(ii) = dispMean;
            end

            % create deltaT array, using a for loop or vectorization
            deltaT = 1:NT;
        end
        
        function [ISF,deltaT] = cal_ISF(obj)            
            if isempty(obj.x_comp)
                obj.cal_c_pos();
            end
            % create MSD array (y-axis of MSD plot)
            ISF = zeros(round(9 * obj.trial.frames/10),1);
            NT = length(ISF);
            q = pi / sqrt(obj.trial.lengthscale(end)*obj.trial.lengthscale(end-1)/(3.14 * obj.trial.Ncell));
            logindex = unique(round(logspace(0, log10(NT),100)));
            % loop over the different possible time windows, calculate MSD for each
            % time window size
            for ii = logindex
                % calculate x displacements, separated by ii indices
                dx = obj.x_comp(obj.trial.offset+ii:end,:) - obj.x_comp(obj.trial.offset:end-ii,:);
                % calculate y displacements similarly
                dy = obj.y_comp(obj.trial.offset+ii:end,:) - obj.y_comp(obj.trial.offset:end-ii,:);
                count = 0;
                isf = 0;
                for th = 0: 0.1: 2*3.141
                % take mean over all displacements
                    isf = isf + mean(real(exp(sqrt(-1) * (cos(th)*q*dx + sin(th)*q*dy))),'all');
                    count = count + 1;
                end
                % store in MSD array
                ISF(ii) = isf/count;
            end
            % create deltaT array, using a for loop or vectorization
            deltaT = 1:NT;
        end
        
        function [xcomp,ycomp] = help_cal_c_pos(obj, xpos_at_frame, ypos_at_frame)
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
        
        function cal_c_pos(obj)
            obj.x_comp=zeros(obj.trial.frames,obj.trial.Ncell);
            obj.y_comp=zeros(obj.trial.frames,obj.trial.Ncell);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.N * ( i - 1 );
                end_point = obj.trial.N * i;
                [xcomp_t,ycomp_t] = obj.help_cal_c_pos(obj.trial.fileReader.coordinate(start_point:end_point,1),obj.trial.fileReader.coordinate(start_point:end_point,2));
                obj.x_comp(i,:)= xcomp_t;
                obj.y_comp(i,:)= ycomp_t;
            end
        end
        
        function cal_c_x_v(obj)
            obj.x_v = zeros(obj.trial.frames,obj.trial.Ncell);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.Ncell * ( i - 1 );
                end_point = obj.trial.Ncell * i;
                obj.x_v(i,:)= obj.trial.fileReader.vel(start_point:end_point,1);
            end
        end
        
        function [tran, rota] = cal_trans_rotat(obj)
            tran = zeros(obj.trial.frames,1);
            rota = zeros(obj.trial.frames,1);
            for i = 1 : obj.trial.frames
                start_point = 1 + obj.trial.Ncell * ( i - 1 );
                end_point = obj.trial.Ncell * i;
                total_K = obj.trial.fileReader.vel(start_point:end_point,4);
                tran_K_fraction = obj.trial.fileReader.vel(start_point:end_point,3);
                tran_K = tran_K_fraction .* total_K;
                tran(i) = mean(tran_K, 'all');
                rota(i) = mean(total_K - tran_K, 'all');
            end
        end
    end
end

