classdef Trial_Disk < Trial_DPM
    %TRIAL_DISK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods
        
        function plotConfig(obj, num, start_point, end_point)
            obj.plot_disk_2d_c(num,[obj.lengthscale(end-1),obj.lengthscale(end)],...
                obj.fileReader.cal_A,obj.fileReader.coordinate(start_point:end_point,1),obj.fileReader.coordinate(start_point:end_point,2))
        end
        
        function plot_disk_2d_c(obj, fig, L, r, x_f, y_f)
            % determine colors
            c       = zeros(obj.Ncell,3);
            c0      = [0 0.2 0.95];
            rmin    = min(obj.lengthscale(1:end-2));
            rs      = obj.lengthscale(1:end-2)./rmin;        
            for n = 1:obj.Ncell
                c(n,:) = rs(n).*c0;
            end
            cmax    = max(max(c));
            c       = c./cmax; 
            c = jet(obj.Ncell);
            c(1,:) = [0,0,0];

            figure(fig), clf, hold on, box on;
            axis('equal');
            axis([0 L(1) 0 L(2)]);

            xcomp = zeros(1,obj.Ncell);
            ycomp = zeros(1,obj.Ncell);
            for ci = 1:obj.Ncell
                index = sum(obj.lengthscale(1:ci),'all');
                start_point_last = 1 + index - obj.lengthscale(ci);
                end_point_last = index;
                cx_tmp = mean(x_f(start_point_last:end_point_last),'all');
                cy_tmp = mean(y_f(start_point_last:end_point_last),'all');
                xcomp(ci) = mod(cx_tmp,obj.lengthscale(end-1));
                ycomp(ci) = mod(cy_tmp,obj.lengthscale(end));
            end

            if (obj.periodic == 1)
                xcomp = mod(xcomp,L(1));
                ycomp = mod(ycomp,L(2));
            end

            for ci = 1:obj.Ncell
                    rectangle('Position',[xcomp(ci)-r(ci), ycomp(ci)-r(ci), 2*r(ci), 2*r(ci)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                    if (obj.periodic == 1)
                        if (xcomp(ci)+r(ci))>L(1)
                            rectangle('Position',[xcomp(ci)-r(ci)-L(1), ycomp(ci)-r(ci), 2*r(ci), 2*r(ci)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end

                        if (xcomp(ci)-r(ci))<0
                            rectangle('Position',[xcomp(ci)-r(ci)+L(1), ycomp(ci)-r(ci), 2*r(ci), 2*r(ci)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end

                        if (ycomp(ci)+r(ci))>L(2)
                            rectangle('Position',[xcomp(ci)-r(ci), ycomp(ci)-r(ci)-L(1), 2*r(ci), 2*r(ci)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end

                        if (ycomp(ci)-r(ci))<0
                            rectangle('Position',[xcomp(ci)-r(ci), ycomp(ci)-r(ci)+L(1), 2*r(ci), 2*r(ci)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end        

                        if ((xcomp(ci)+r(ci))>L(1) || (xcomp(ci)-r(ci))<0) && ((ycomp(ci)+r(ci))>L(2) || (ycomp(ci)-r(ci))<0)
                            x1 = mod((xcomp(ci) + L(1)),L(1));
                            y1 = mod((ycomp(ci) + L(2)),L(2));
                            rectangle('Position',[x1-r(ci), y1-r(ci), 2*r(ci), 2*r(ci)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end
                    end
            end
            drawnow;
            %hold off;
        end
    end
end

