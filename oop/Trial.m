classdef Trial < handle
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
    end
    
    methods
        function obj = Trial(t_index_i,t_index_j,basefolder,fileReader)
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
        
        function plotConfig(obj, num, start_point, end_point)
            obj.plot_particles_2d_c(num,[obj.lengthscale(end-1),obj.lengthscale(end)],...
                obj.fileReader.coordinate(start_point:end_point,3)/2,obj.fileReader.coordinate(start_point:end_point,1),obj.fileReader.coordinate(start_point:end_point,2))
        end
        
        function createCalculator(obj)
            obj.calculator = Calculator(obj);
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

