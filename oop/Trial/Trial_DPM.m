classdef Trial_DPM < handle
    %TRIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        basefolder
        t_index_i
        t_index_j
        folder
        fileReader
        periodic = 1
        %coordinate
        N
        frames
        Ncell
        lengthscale
        %vel
        calculator
        velocity
        Dr
        kb
        kl
        calA0
        MSD
        deltaT
        logindex
    end
    
    methods
        function obj = Trial_DPM(t_index_i,t_index_j,basefolder,fileReader)
            %TRIAL Construct an instance of this class
            %   Detailed explanation goes here
            obj.t_index_i = t_index_i;
            obj.t_index_j = t_index_j;
            obj.basefolder = basefolder;
            obj.fileReader = fileReader;
            obj.fileReader.setFolder(basefolder, t_index_i, t_index_j);
        end
        
        function getMDinfo(obj)
            obj.lengthscale = obj.fileReader.lengthscale;
            obj.N=sum(obj.lengthscale(1:end-2),'all');
            obj.frames= size(obj.fileReader.coordinate,1)/obj.N ;
            obj.Ncell = size(obj.lengthscale,1)-2;
        end
        
        function readInitial(obj)
            obj.fileReader.readInitial();
            %obj.coordinate = obj.fileReader.coordinate;
            obj.getMDinfo();
        end
        
        function readV0(obj)
            obj.fileReader.readV0();
            obj.velocity = obj.fileReader.v0(1);
            obj.Dr = obj.fileReader.v0(2);
            obj.kb = obj.fileReader.v0(3);
            obj.kl = obj.fileReader.v0(4);
            obj.calA0 = obj.fileReader.v0(5);
        end
        
        function plotInitial(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.readInitial();
            obj.plotLastFrame(1);
        end
        
        function readMDdata(obj)
            obj.fileReader.readMDdata();
            %obj.coordinate = obj.fileReader.coordinate;
            %obj.vel = obj.fileReader.vel;
            obj.getMDinfo();
        end
        
        function plotLastFrame(obj, num)
            i = floor(obj.frames);
            start_point = 1 + obj.N * ( i - 1 );
            end_point = obj.N * i;
            obj.plotConfig(num, start_point, end_point);
        end
        
        function plotVelDistribution(obj)
            vel_s = sqrt(obj.fileReader.vel(:,1).^2 + obj.fileReader.vel(:,2).^2)/obj.fileReader.lengthscale(end);
            %vel_s = obj.fileReader.vel(:,2)/obj.fileReader.lengthscale(end);
            figure(8), hold on, box on;
            set(gcf,'color','w');
            h = histogram(vel_s,'BinWidth',5e-5, 'Normalization', 'probability');
        end
        
        function plotCalADistribution(obj)
            figure(9), hold on, box on;
            set(gcf,'color','w');
            h = histogram(obj.fileReader.cal_A,'BinWidth',0.01, 'Normalization', 'probability');
        end
        
        function plotRotationVsTranslaion(obj)
            [tran, rota] = obj.calculator.cal_trans_rotat();
            figure(6), hold on, box on;
            set(gcf,'color','w');
            scatter(1:obj.frames, tran);
            scatter(1:obj.frames, rota);
            yline(mean(tran,'all'));
            yline(mean(rota,'all'));
            hold off;
        end
        
        function plotConfig(obj, num, start_point, end_point)
            obj.plot_particles_2d_c(num,[obj.lengthscale(end-1),obj.lengthscale(end)],...
                obj.fileReader.coordinate(start_point:end_point,3)/2,obj.fileReader.coordinate(start_point:end_point,1),obj.fileReader.coordinate(start_point:end_point,2))
        end
        
        function createCalculator(obj)
            obj.calculator = Calculator(obj);
        end
        
        function rescaleMSD(obj)
            obj.MSD = obj.MSD/(obj.fileReader.lengthscale(end)^2);
            time_scale = (1/5000) * (100000/0.005);
            obj.deltaT = obj.deltaT * time_scale;
        end
        
        function cal_msd(obj)
            [obj.MSD, obj.deltaT] = obj.calculator.cal_msd();
            obj.logindex = unique(round(logspace(0, log10(obj.deltaT(end)),1000)));
            obj.rescaleMSD();
        end
        
        function plotMSD(obj)
            figure((obj.t_index_j+1)*10+2)
            set(gcf,'color','w');
            hold on, box on;
            % plot curve, add units to axes, etc
            %plot(deltaT(logindex), MSD(logindex),'color','red','linewidth',3);
            
            plot(obj.deltaT(obj.logindex), obj.MSD(obj.logindex),'linewidth',3);
            xlabel('time');ylabel('MSD');
            length_t = length(obj.deltaT);
            half_log = round(size(obj.logindex,2)/3);
            if half_log == 0
                return
            end
            [P,S] = polyfit(log10(obj.deltaT(obj.logindex(half_log:end))), log10(obj.MSD(obj.logindex(half_log:end))'), 1);
            uncertainty = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
            %P = polyfit(log10(deltaT(1: end)), log10(MSD(1: end))', 1);
            yfit = P(1)*log10(obj.deltaT(obj.logindex))+P(2);
            %plot(obj.deltaT(obj.logindex),10.^(yfit),'r-.');
            theString = sprintf('slope = %.3f ', P(1));
            %text(10^5, 0.01, theString, 'FontSize', 20);
            ax = gca;
            %ax.FontSize = 22;
            ax.XScale = "log";
            ax.YScale = "log";
            
            disp(["difusion: ", (10^P(2))/(4 * 0.005), " +- ", (10^P(2) * log(10) * uncertainty(2))/(4*0.005)])
            disp(["slope: ", P(1), "+- ", uncertainty(1)])
            disp(["intercept: ", P(2), "+- ", uncertainty(2)])
        end
        
        function showVideo(obj, tFrame)
            for i = 1 : ceil(obj.frames/tFrame) : obj.frames
                start_point = 1 + obj.N * ( i - 1 );
                end_point = obj.N * i;
                obj.plotConfig(2, start_point, end_point);
            end
        end
        
        function saveVideo(obj, tFrame)
            loc = strfind(obj.basefolder,'/');
            name = extractAfter(obj.basefolder,loc(2));
            name = extractBefore(name,'/');
            vobj = VideoWriter("video/" + name + obj.t_index_i + obj.t_index_j + ".avi");
            vobj.FrameRate = 5;
            open(vobj);
            for i = 100 : ceil(obj.frames/tFrame) : obj.frames
                start_point = 1 + obj.N * ( i - 1 );
                end_point = obj.N * i;
                obj.plotConfig(2, start_point, end_point);
                set(gcf,'color','w');
                frame = getframe(gcf) ;
                writeVideo(vobj, frame);
            end
            close(vobj);
        end
        
        function plot_particles_2d_c(obj, fig, L, r_f, x_f, y_f)
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

            if (obj.periodic == 1)
                axis([0 L(1) 0 L(2)]);
                x_f = mod(x_f,L(1));
                y_f = mod(y_f,L(2));
            else
                axis([0 L(1)*1.1 0 L(2)]);
            end

            for ci = 1:obj.Ncell
                index = sum(obj.lengthscale(1:ci),'all');
                start_point_last = 1 + index - obj.lengthscale(ci);
                end_point_last = index;
                x = x_f(start_point_last:end_point_last);
                y = y_f(start_point_last:end_point_last);
                r = r_f(start_point_last:end_point_last);
                for n = 1:obj.lengthscale(ci)

                    rectangle('Position',[x(n)-r(n), y(n)-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                    if (obj.periodic == 1)
                        if (x(n)+r(n))>L(1)
                            rectangle('Position',[x(n)-r(n)-L(1), y(n)-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end

                        if (x(n)-r(n))<0
                            rectangle('Position',[x(n)-r(n)+L(1), y(n)-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end

                        if (y(n)+r(n))>L(2)
                            rectangle('Position',[x(n)-r(n), y(n)-r(n)-L(1), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end

                        if (y(n)-r(n))<0
                            rectangle('Position',[x(n)-r(n), y(n)-r(n)+L(1), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end        

                        if ((x(n)+r(n))>L(1) || (x(n)-r(n))<0) && ((y(n)+r(n))>L(2) || (y(n)-r(n))<0)
                            x1 = mod((x(n) + L(1)),L(1));
                            y1 = mod((y(n) + L(2)),L(2));
                            rectangle('Position',[x1-r(n), y1-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(ci,:),'facecolor',c(ci,:));
                        end
                    end

                end
            end
            drawnow;
            %hold off;
        end
    end
end

