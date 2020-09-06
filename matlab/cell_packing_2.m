clc
clear
close all

global T1_count;
global T1_cells;
global T1_index;

mode = 1;

%basefolder = "D:\project\cells1_N\";
%basefolder = "D:\project\cells38\";
basefolder = "~/project/cells54/";
%basefolder = "~/project/test/";
%basefolder = "~/scratch60/cells47/";
%basefolder = "~/project/cells48/";
%basefolder = "C:\Users\Yuxuan Cheng\source\repos\cells\forked-cells\forked-cells\";

all_mean_cal_A = [];
order_per = [];
ifjammed = [];
mean_p = [];
var_p = [];
msd = {};
ISF_en = {};
tao_en = {};
T1_cells_list = {};
v0 = [];
v0_en = {};

time_scale = (1/5000) * (100000/0.005);

for t_index_i =2:9
    %close all
    for t_index_j = 0:9
    try
    
    if(mode==1)
        folder = basefolder + int2str(t_index_i) + "/";
        v0_file = folder + "v0.txt";
        tmp = csvread(v0_file);
        v0_en{t_index_i+1,t_index_j+1} = tmp(t_index_j+1,:);
        v0 = [v0;csvread(v0_file)];
    elseif(mode==2)
        folder = basefolder + int2str(t_index_i) + "_" + int2str(t_index_j) + "/";
        v0_file = folder + "v0.txt";
        tmp = csvread(v0_file);
        v0_en{t_index_i+1,t_index_j+1} = tmp;
        v0 = [v0;csvread(v0_file)];
    elseif(mode==0)
        folder = basefolder;
        v0_file = folder + "v0.txt";
        v0 = csvread(v0_file);
    end
    v0 = unique(v0,'rows');
    %extend = "_jammed_" + int2str(t_index_i) + int2str(t_index_j) +".txt";
    extend = "_jammed_" + int2str(t_index_i) +".txt";
    %extend = ".txt";
    
    T1_cells = [];
    T1_index = [];
    
    coordinate_file = folder + "jam" + extend;
    length_file = folder + "length" + extend;
    cal_A_file = folder + "calA" + extend;
    contact_file = folder + "contact" + extend;
    v_file = folder + "v" + extend;  
        
    packing_contact = dlmread(contact_file);
    packing_contact = packing_contact(end,:);
    coordinate = csvread(coordinate_file);
    %length_file = folder + "length" + extend;
    lengthscale = csvread(length_file);
    cal_A = csvread(cal_A_file);
    N=sum(lengthscale(1:end-2),'all');
    frames= size(coordinate,1)/N ;
    Ncell = size(lengthscale,1)-2;

%     for i = 1 :frames
%         start_point = 1 + N * ( i - 1 );
%         end_point = N * i;
%         plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%             coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%     %         frame = getframe(gcf) ;
%     %         writeVideo(vobj, frame);
%     end

    i = frames;
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d(1,[lengthscale(end-1),lengthscale(end)],...
    coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))

    xpos_at_frame = coordinate(start_point:end_point,1);
    ypos_at_frame = coordinate(start_point:end_point,2);
    [vAll,cAll,xcomp_j,ycomp_j] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, 3);
    jammed_voronoi_con_net = voronoi_contact(cAll, Ncell); 
        
        
    extend1 = "_" + int2str(t_index_i) + int2str(t_index_j) +".txt";

    coordinate_file = folder + "jam" + extend1;
    coordinate = csvread(coordinate_file);
    %length_file = folder + "length" + extend;
    length_file = folder + "length" + extend;
    lengthscale = csvread(length_file);
    cal_A_file = folder + "calA" + extend1;
    cal_A = csvread(cal_A_file);
    try
        contact_file = folder + "contact" + extend1;
        contact = dlmread(contact_file);
    catch
        disp('no contact file')
    end
    v_file = folder + "v" + extend1;
    vel = csvread(v_file);
    
    catch
        continue
    end
    
%     contact = contact(end,:);
%     N=sum(lengthscale(1:end-2),'all');
%     frames= size(coordinate,1)/N ;
%     difference = sum(contact ~= packing_contact,'all');
%     if difference > 0
%         disp(t_index)
%         disp(difference)
%     end
    
    N=sum(lengthscale(1:end-2),'all');    
        
%     frames= floor((size(coordinate,1)/N)/2);
%     coordinate = coordinate(frames*N + 1:end,:);
    
    frames= size(coordinate,1)/N ;
    %frames = 3700;
    Ncell = size(lengthscale,1)-2;


%     %vobj = VideoWriter('test.mp4','MPEG-4');
%     vobj = VideoWriter('test1.avi');
%     vobj.FrameRate = 3;
%     open(vobj);
%         for i = 1 :  round(frames/50):frames
%             start_point = 1 + N * ( i - 1 );
%             end_point = N * i;
%             plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%                 coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%             frame = getframe(gcf) ;
%             writeVideo(vobj, frame);
%         end
% 
%     close(vobj);


%     for i = 1 :  round(frames/20):frames
%         start_point = 1 + N * ( i - 1 );
%         end_point = N * i;
%         plot_particles_2d(2,[lengthscale(end-1),lengthscale(end)],...
%             coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
%     end
% % %     
%     T1_count = zeros(1, frames);
%     for i = 1 : frames
%         start_point = 1 + N * ( i - 1 );
%         end_point = N * i;
%         xpos_at_frame = coordinate(start_point:end_point,1);
%         ypos_at_frame = coordinate(start_point:end_point,2);
%         T1_count(1,i) = T1_swap(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, (t_index_j+1)*10+1);
%     end
%     

%     T1_count = zeros(1, frames);
%     for i = 1 : frames-1
%         start_point = 1 + N * ( i - 1 );
%         end_point = N * i;
%         xpos_at_frame = coordinate(start_point:end_point,1);
%         ypos_at_frame = coordinate(start_point:end_point,2);
%         start_point = 1 + N * ( i );
%         end_point = N * (i+1);
%         xpos_at_next_frame = coordinate(start_point:end_point,1);
%         ypos_at_next_frame = coordinate(start_point:end_point,2);
%         %T1_count(1,i+1) = T1_swap_1(xpos_at_frame, ypos_at_frame, ...
%         %xpos_at_next_frame, ypos_at_next_frame, Ncell,lengthscale, ...
%         %(t_index_j+1)*10+1);
%         %% 1 for rate, 2 for cumulated T1, 3 for cell fraction
%         T1_count(1,i+1) = T1_swap_2(xpos_at_frame, ypos_at_frame, ...
%         xpos_at_next_frame, ypos_at_next_frame, Ncell,lengthscale, i, ...
%         (t_index_j+1)*10+1,3);
%     end
%     
%     T1_cells_list{t_index_i+1,t_index_j+1} = T1_count;
%     
%     figure((t_index_j+1)*10+5), clf, hold on, box on;
%     n = 100; % average every n values
%     %T1_count = arrayfun(@(i) mean(T1_count(i:i+n-1)),1:n:length(T1_count)-n+1)'; % the averaged vector
%     plot((1:length(T1_count))* time_scale,T1_count);
%     xlabel('time');ylabel('T1');
%     ax = gca;
%     ax.XScale = "log";
%     hold off;
%     %sum(T1_count)/(frames *(1/5000) * (100000/0.005))
%     figure((t_index_j+1)*10+6)
%     test = find(T1_count>0) * (1/5000) * (100000/0.005);
%     [count,edges] = histcounts(log10(test));
%     histogram(test,10.^edges,'Normalization','countdensity')
%     xlabel('time');ylabel('T1 rate');
%     set(gca, 'xscale','log')
%     
%     i = frames;
%     
%         vobj = VideoWriter('test.mp4','MPEG-4');
%         vobj.FrameRate = 10;
%         open(vobj);
%             for i = 1 :  round(frames/100):frames
%                 start_point = 1 + N * ( i - 1 );
%                 end_point = N * i;
%                 xpos_at_frame = coordinate(start_point:end_point,1);
%                 ypos_at_frame = coordinate(start_point:end_point,2);
%                 [vAll,cAll,xcomp,ycomp] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, (t_index_j+1)*10+1);
%                 frame = getframe(gcf) ;
%                 writeVideo(vobj, frame);
%             end
%         
%         close(vobj);

    i = floor(frames);
    start_point = 1 + N * ( i - 1 );
    end_point = N * i;
    plot_particles_2d((t_index_j+1)*10,[lengthscale(end-1),lengthscale(end)],...
        coordinate(start_point:end_point,3)/2,coordinate(start_point:end_point,1),coordinate(start_point:end_point,2))
    xpos_at_frame = coordinate(start_point:end_point,1);
    ypos_at_frame = coordinate(start_point:end_point,2);
    [vAll,cAll,xcomp,ycomp] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, (t_index_j+1)*10+1);
    voronoi_con_net = voronoi_contact(cAll, Ncell);
    
%     [structure, wave] = cal_structure_factor(N, coordinate, frames, Ncell, lengthscale);   
%     figure(5)
%     plot(wave, structure);
%     xlabel('q');ylabel('S');
    
    %is_equal = isequal(voronoi_con_net,jammed_voronoi_con_net);
    %is_equal = compare_contact(voronoi_con_net,jammed_voronoi_con_net,Ncell);
    %is_equal = compare_contact1(voronoi_con_net,jammed_voronoi_con_net,xcomp_j,ycomp_j,lengthscale,Ncell);
    
    [MSD,deltaT] = cal_msd(N, coordinate, frames, Ncell, lengthscale);
    %[MSD,deltaT] = cal_msd_vertex(N, coordinate, frames);
    deltaT = deltaT * time_scale;
    
    % open figure window
    figure((t_index_j+1)*10+2), clf, hold on, box on;
    % plot curve, add units to axes, etc
    plot(deltaT, MSD,'color','red','linewidth',3);
    xlabel('time');ylabel('MSD');
    length_t = length(deltaT);
    P = polyfit(log10(deltaT(round(3*length_t/6): end)), log10(MSD(round(3*length_t/6): end))', 1);
    %P = polyfit(log10(deltaT(1: end)), log10(MSD(1: end))', 1);
    yfit = P(1)*log10(deltaT)+P(2);
    plot(deltaT,10.^(yfit),'r-.');
    theString = sprintf('slope = %.3f ', P(1));
    text(10^5, 0.01, theString, 'FontSize', 20);
    ax = gca;
    %ax.FontSize = 22;
    ax.XScale = "log";
    ax.YScale = "log";
    
    msd{t_index_i+1,t_index_j+1} = MSD;
    cri = max(MSD,[],'all')/(lengthscale(end)*lengthscale(end-1)/(3.14*Ncell));
    is_equal = (max(MSD,[],'all')<(lengthscale(end)*lengthscale(end-1)/(3.14*Ncell)));
    
    ifjammed = [ifjammed, cri];
    %ifjammed = [ifjammed, is_equal];
    
    if is_equal == 0
        disp({'voronoi_con_net is different for ',t_index_i,t_index_j})
    end
    
    [ISF,deltaT1] = cal_ISF(N, coordinate, frames, Ncell, lengthscale);
    logindex = unique(round(logspace(0, log10(deltaT1(end)),100)));
    deltaT1 = deltaT1 * time_scale;  
    
%     fitfun = fittype( @(C, tao, b, x) C*exp(-(x/tao).^b));
%     fitted = fit( (logindex)', ISF(logindex), fitfun, 'StartPoint', [1,500,1]);
%     
    fitfun = fittype( @(tao, b, x) exp(-(x/tao).^b));
    fitted = fit( (logindex)', abs(ISF(logindex)), fitfun, 'StartPoint', [500,1]);
    
    ISF_en{t_index_i+1,t_index_j+1} = ISF;
    tao_en{t_index_i+1,t_index_j+1} = fitted.tao * time_scale;
    % open figure window
    figure(13), hold on, box on;
    if t_index_j == 0
        clf
    end
    %figure((t_index_j+1)*10+3), clf, hold on, box on;
    % plot curve, add units to axes, etc
    %plot(deltaT1(logindex), ISF(logindex),'color','red','linewidth',3);
    plot(fitted,(logindex), ISF(logindex));
    xlabel('time');ylabel('ISF');
    length_t = length(deltaT);
    ax = gca;
    ax.XScale = "log";
    ax.YScale = "Linear";
    
    
    v_d_index = t_index_i * 10 + t_index_j + 1;
    
%     order = cal_order1(vel)/v0(v_d_index,1);
%     order_per = [order_per, order];
    
%     
%     n_particle = density_fluctuation(N, coordinate, frames, Ncell, lengthscale);
%     N_var = var(n_particle);
%     mean_N = mean(n_particle);
%     disp(N_var/sqrt(mean_N))
%     mean_p = [mean_p, mean_N];
%     var_p = [var_p, N_var];
    
%     [std_n, mean_n] = giant_number(N, coordinate, frames, Ncell, lengthscale);
%     
%     figure((t_index_j+1)*10+4);hold on;
% %     mean_n = mean_n(3:end);std_n = std_n(3:end);
%     % loglog(mean_p,std_p,'-s')
%     scatter(mean_n,std_n,25,'b','*') 
%     P = polyfit(log10(mean_n),log10(std_n),1);
%     yfit = P(1)*log10(mean_n)+P(2);
%     plot(mean_n,10.^yfit,'r-.');
%     theString = sprintf('slope = %.3f ', P(1));
%     text(5, 5, theString, 'FontSize', 24);
%     xlabel("Mean(N)");
%     ylabel("Std(N)"); 
%     ax = gca;
%     ax.XScale = "log";
%     ax.YScale = "log";
    
    
    cal_A = reshape(cal_A,Ncell,[]);
    mean_cal_A = mean(cal_A,'all');
    
    all_mean_cal_A = [all_mean_cal_A, mean_cal_A];
    
    disp({'cell index',t_index_i,t_index_j})
%    [breaked, new] = compare_contact2(voronoi_con_net,jammed_voronoi_con_net,Ncell)
%     sum((round(mean(contact,1))-packing_contact)==-1)
%     sum((round(mean(contact,1))-packing_contact)==1)
%     sum((contact(end,:)-packing_contact)==-1)
%     sum((contact(end,:)-packing_contact)==1)
    
    end
end

% order_per = reshape(order_per, 10, []);

% figure(5);
% heatmap(order_per)
% xlabel("v0");
% ylabel("noise");
% title("order parameter");

% figure(7), clf, hold on, box on;
% plot(msd{9,3},'color','red','linewidth',3);
% plot(msd{9,6},'color','green','linewidth',3);
% plot(msd{9,7},'color','black','linewidth',3);
% plot(msd{9,10},'color','blue','linewidth',3);
% legend({'Liquid like','Solid like','fit'});
% xlabel('time');ylabel('MSD');
% ax = gca;
% ax.FontSize = 22;
% ax.XScale = "log";
% ax.YScale = "log";

figure(7), clf, hold on, box on;
% plot curve, add units to axes, etc
%deltaT = deltaT * 1/0.005;
plot(deltaT, msd{9,3},'color','red','linewidth',3);
plot(deltaT, msd{9,6},'color','green','linewidth',3);
plot(deltaT, msd{9,7},'color','black','linewidth',3);
plot(deltaT, msd{9,10},'color','blue','linewidth',3);
%plot(deltaT, msd{3,10},'color','blue','linewidth',3);
%P = polyfit(log10(deltaT(round(5*frames/6): end)), log10(MSD(round(5*frames/6): end))', 1);
MSD1=msd{9,10};
P = polyfit(log10(deltaT(round(3*length_t/6): end)), log10(MSD1(round(3*length_t/6): end))', 1);
yfit = P(1)*log10(deltaT)+P(2);
plot(deltaT,10.^(yfit),'r-.');
legend({'Liquid like','Solid like','fit'});
xlabel('time');ylabel('MSD');
theString = sprintf('slope = %.3f ', P(1));
%text(5, 1, theString, 'FontSize', 24);
ax = gca;
ax.FontSize = 22;
ax.XScale = "log";
ax.YScale = "log";

figure(200), clf, hold on, box on;
plot(deltaT(logindex), ISF_en{7,1}(logindex),'color','red','linewidth',3);
plot(deltaT(logindex), ISF_en{7,3}(logindex),'color','green','linewidth',3);
plot(deltaT(logindex), ISF_en{7,6}(logindex),'color','black','linewidth',3);
plot(deltaT(logindex), ISF_en{7,10}(logindex),'color','blue','linewidth',3);
xlabel('time');ylabel('ISF');
ax = gca;
ax.FontSize = 22;
ax.XScale = "log";
% ax.YScale = "log";

figure(201), clf, hold on, box on;
% plot curve, add units to axes, etc
deltaT1 = (1:length(T1_count))* time_scale;
% for i = 1 : 10
%     plot(deltaT1, T1_cells_list{6,i});
% end
plot(deltaT1, T1_cells_list{6,1},'color','red');
plot(deltaT1, T1_cells_list{6,5},'color','green');
plot(deltaT1, T1_cells_list{6,6},'color','black');
plot(deltaT1, T1_cells_list{6,10},'color','blue');
xlabel('time');ylabel('T1');
ax = gca;
ax.XScale = "log";
ax.YScale = "log";

% all_mean_cal_A = reshape(all_mean_cal_A, [], 10);
% figure(8);
% heatmap(flip(all_mean_cal_A,1))
% ax = gca;
% ax.XData = unique(v0(:,end-1));
% ax.YData = flip(unique(v0(:,1)));
% xlabel("calA0");
% ylabel("v0");
% title("calA");

% all_mean_cal_A = reshape(all_mean_cal_A, [], 10);
% figure(8);
% heatmap(flip(all_mean_cal_A,1))
% ax = gca;
% ax.XData = unique(v0(:,3));
% ax.YData = flip(unique(v0(:,1)));
% xlabel("calA0");
% ylabel("v0");
% title("calA");
% 
% ifjammed = reshape(ifjammed, [], 10);
% figure(6);
% %heatmap(flip(ifjammed,1),'CellLabelColor','none')
% heatmap(flip(ifjammed,1))
% ax = gca;
% ax.XData = unique(v0(:,end-1));
% ax.YData = flip(unique(v0(:,1)));
% xlabel("calA0");
% ylabel("v0");
% title("Phase Diagram");

ifjammed = reshape(ifjammed, [], 10);
figure(6);
%heatmap(flip(ifjammed,1),'CellLabelColor','none')
heatmap(flip(ifjammed,1))
ax = gca;
ax.XData = unique(v0(:,3));
ax.YData = flip(unique(v0(:,1)));
xlabel("kb");
ylabel("v0");
title("Phase Diagram");

all_mean_cal_A = reshape(all_mean_cal_A, [], 10);
figure(8);
%heatmap(flip(ifjammed,1),'CellLabelColor','none')
heatmap(flip(all_mean_cal_A,1))
ax = gca;
ax.XData = unique(v0(:,3));
ax.YData = flip(unique(v0(:,1)));
xlabel("kb");
ylabel("v0");
title("Phase Diagram");

figure(9);
all_mean_cal_A = reshape(all_mean_cal_A, [], 10);
re = all_mean_cal_A';
scatter(1./(v0(:,1).^2),[tao_en{:}].*v0(:,1)',30,re(:),'filled');
cb = colorbar();
ax = gca;
ax.XScale = "log";
ax.YScale = "log";
xlabel('1/v^2');ylabel('tao * v');

% figure(9); hold on
% re = all_mean_cal_A';
% for i = 1:10
%     plot(1./(v0(i:10:end,1).^2),[tao_en{i,:}].*v0(i:10:end,1)','color','blue');
% end
% scatter(1./(v0(:,1).^2),[tao_en{:}].*v0(:,1)',30,re(:),'filled');
% cb = colorbar();
% %xlim([-inf, 10^4]);
% ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
% xlabel('1/v^2');ylabel('tao * v');

figure(9); hold on
v0_s = [];
re = all_mean_cal_A';
for i = 1:10
    v0_temp = cell2mat(v0_en(i,:)');
    v0_temp = sort(v0_temp(:,1));
    v0_s = [v0_s; v0_temp'];
    plot(1./(v0_temp.^2),[tao_en{i,:}].*v0_temp');
end
% scatter(1./(v0_s(:).^2),[tao_en{:}].*v0_s(:)',30,re(:),'filled');
% cb = colorbar();
% %xlim([-inf, 10^4]);
ax = gca;
ax.XScale = "log";
ax.YScale = "log";
xlabel('1/v^2');ylabel('tao * v');

figure(9); hold on
v0_s = [];
re = all_mean_cal_A';
for i = 3:10
    v0_temp = cell2mat(v0_en(i,2:end)');
    v0_temp = sort(v0_temp(:,1));
    v0_s = [v0_s; v0_temp'];
    plot(1./(v0_temp.^2),[tao_en{i,2:end}].*v0_temp');
end
% scatter(1./(v0_s(:).^2),[tao_en{:}].*v0_s(:)',30,re(:),'filled');
% cb = colorbar();
% %xlim([-inf, 10^4]);
ax = gca;
ax.XScale = "log";
ax.YScale = "log";
xlabel('1/v^2');ylabel('tao * v');

figure(9); hold on
kb_0 = 0.00037;
u = 10;
sigma = 0.005;
v0_s = [];
re = all_mean_cal_A';
for i = 3:10
    v0_temp = cell2mat(v0_en(i,2:end)');
    kb_temp = v0_temp(1,3);
    v0_temp = sort(v0_temp(:,1));
    v0_s = [v0_s; v0_temp'];
    x = abs(kb_temp - kb_0);
    plot((x^(2/u))./(v0_temp.^2),x^sigma * log10([tao_en{i,2:end}].*v0_temp'));
end
% scatter(1./(v0_s(:).^2),[tao_en{:}].*v0_s(:)',30,re(:),'filled');
% cb = colorbar();
% %xlim([-inf, 10^4]);
ax = gca;
ax.XScale = "log";
ax.YScale = "log";
xlabel('1/v^2');ylabel('tao * v');

% all_mean_cal_A = reshape(all_mean_cal_A, 10, []);
% figure(8);
% heatmap(flip(all_mean_cal_A,1))
% ax = gca;
% ax.XData = unique(v0(:,end));
% ax.YData = flip(unique(v0(:,1)));
% xlabel("phi");
% ylabel("v0");
% title("calA");
% 
% ifjammed = reshape(ifjammed, 10, []);
% figure(6);
% heatmap(flip(ifjammed,1),'CellLabelColor','none')
% ax = gca;
% ax.XData = unique(v0(:,end));
% ax.YData = flip(unique(v0(:,1)));
% xlabel("phi");
% ylabel("v0");
% title("Phase Diagram");


% figure(6);hold on;
% fliped_jam = flip(ifjammed,1).*(-1) + 1;
% imagesc(fliped_jam)
% [X, Y]=meshgrid(1:10,1:10);
% string = mat2cell(num2str(reshape(flip(all_mean_cal_A),[100,1])),ones(10*10,1));
% text(Y(:)-.5,X(:)+.25,string,'HorizontalAlignment','left')


% figure(7);
% plot(order_per(:,1))
% ylabel("order");
% xlabel("noise");

% figure(7);hold on;
% % loglog(mean_p,std_p,'-s')
% scatter(log10(mean_p),log10(var_p),25,'b','*') 
% P = polyfit(log10(mean_p),log10(var_p),1);
% yfit = P(1)*log10(mean_p)+P(2);
% plot(log10(mean_p),yfit,'r-.');
% theString = sprintf('slope = %.3f ', P(1));
% text(1, 0.4, theString, 'FontSize', 24);
% xlabel("log Mean");
% ylabel("log var"); 

% v0_file = folder + "v0.txt";
% v0 = csvread(v0_file);
% figure(3);
% plot(v0,all_mean_cal_A)
% xlabel("v0");
% ylabel("CalA");



function [vAll,cAll,xcomp,ycomp] = draw_voronoi(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, fig)
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX = [];
    voronoiY = [];
    for i = -1:1
        for j = -1:1
            voronoiX = [voronoiX,xcomp + lengthscale(end-1) * i];
            voronoiY = [voronoiY,ycomp + lengthscale(end)* j];
        end
    end
    text_range = (1 + 4 * Ncell) : (5 * Ncell);
    figure(fig); hold on
    clf;
    voronoi(voronoiX,voronoiY)
    plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell))}, text_range');
     text(voronoiX(text_range), voronoiY(text_range), plabels, 'FontWeight', ...
           'bold', 'HorizontalAlignment','center', ...
           'BackgroundColor', 'none')
    xlim([0,lengthscale(end-1)]);
    ylim([0,lengthscale(end)]);
    hold off
    dt = delaunayTriangulation(voronoiX',voronoiY');
    [vAll,cAll] = voronoiDiagram(dt);
end

% function count = T1_swap(xpos_at_frame, ypos_at_frame, Ncell, lengthscale, fig)
% 
%     threshold = (1/100) * sqrt(lengthscale(end)*lengthscale(end-1)/(3.14*Ncell));
%     xcomp = zeros(1,Ncell);
%     ycomp = zeros(1,Ncell);
%     for ci = 1:Ncell
%         index = sum(lengthscale(1:ci),'all');
%         start_point_last = 1 + index - lengthscale(ci);
%         end_point_last = index;
%         cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
%         cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
%         xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
%         ycomp(ci) = mod(cy_tmp,lengthscale(end));
%     end
%     voronoiX = [];
%     voronoiY = [];
%     for i = -1:1
%         for j = -1:1
%             voronoiX = [voronoiX,xcomp + lengthscale(end-1) * i];
%             voronoiY = [voronoiY,ycomp + lengthscale(end)* j];
%         end
%     end
%     text_range = (1 + 4 * Ncell) : (5 * Ncell);
%     [vx,vy] = voronoi(voronoiX,voronoiY);
%     dist = sqrt((vx(1,:)-vx(2,:)).^2 + (vy(1,:)-vy(2,:)).^2);
%     %count = round(sum(dist < threshold,'all')/9);
%     index = find(dist < threshold);
%     count = sum(vx(1,index)>0&vx(1,index)<lengthscale(end-1)&vy(1,index)>0&vy(1,index)<lengthscale(end),'all');
% %     if count > 0
% %         figure(fig); hold on
% %         clf;
% %         voronoi(voronoiX,voronoiY)
% %         plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell))}, text_range');
% %          text(voronoiX(text_range), voronoiY(text_range), plabels, 'FontWeight', ...
% %                'bold', 'HorizontalAlignment','center', ...
% %                'BackgroundColor', 'none')
% %         xlim([0,lengthscale(end-1)]);
% %         ylim([0,lengthscale(end)]);
% %         hold off
% %     end
% end

function count = T1_swap_1(xpos_at_frame, ypos_at_frame, xpos_at_next_frame, ypos_at_next_frame, Ncell, lengthscale, fig)
    count = 0;
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX = [];
    voronoiY = [];
    for i = -1:1
        for j = -1:1
            voronoiX = [voronoiX,xcomp + lengthscale(end-1) * i];
            voronoiY = [voronoiY,ycomp + lengthscale(end)* j];
        end
    end

    dt = delaunayTriangulation(voronoiX',voronoiY');
    [vAll,cAll] = voronoiDiagram(dt);
    
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_next_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_next_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX_n = [];
    voronoiY_n = [];
    for i = -1:1
        for j = -1:1
            voronoiX_n = [voronoiX_n,xcomp + lengthscale(end-1) * i];
            voronoiY_n = [voronoiY_n,ycomp + lengthscale(end)* j];
        end
    end

    dt = delaunayTriangulation(voronoiX_n',voronoiY_n');
    [vAll_next,cAll_next] = voronoiDiagram(dt);
    
    for i = (1 + 4 * Ncell) : (5 * Ncell)
        if length(cAll{i})~=length(cAll_next{i})
            count = count + 1;
            %disp(mod(i-1,Ncell));
        end
    end
    
%      if mod(count,4) > 0
%         text_range = (1 + 4 * Ncell) : (5 * Ncell);
%         figure(fig); hold on
%         clf;
%         voronoi(voronoiX,voronoiY)
%         plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell))}, text_range');
%          text(voronoiX(text_range), voronoiY(text_range), plabels, 'FontWeight', ...
%                'bold', 'HorizontalAlignment','center', ...
%                'BackgroundColor', 'none')
%         xlim([0,lengthscale(end-1)]);
%         ylim([0,lengthscale(end)]);
%         hold off
%         
%         figure(fig+1); hold on
%         clf;
%         voronoi(voronoiX_n,voronoiY_n)
%         plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell))}, text_range');
%          text(voronoiX_n(text_range), voronoiY_n(text_range), plabels, 'FontWeight', ...
%                'bold', 'HorizontalAlignment','center', ...
%                'BackgroundColor', 'none')
%         xlim([0,lengthscale(end-1)]);
%         ylim([0,lengthscale(end)]);
%         hold off
%      end
    count = ceil(count / 4);
end

function count = T1_swap_2(xpos_at_frame, ypos_at_frame, xpos_at_next_frame, ypos_at_next_frame, Ncell, lengthscale, frame_i, fig, mode)
    global T1_cells;
    global T1_index;
    global T1_count;
    count = 0;
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX = [];
    voronoiY = [];
    for i = -1:1
        for j = -1:1
            voronoiX = [voronoiX,xcomp + lengthscale(end-1) * i];
            voronoiY = [voronoiY,ycomp + lengthscale(end)* j];
        end
    end

    dt = delaunayTriangulation(voronoiX',voronoiY');
    [vAll,cAll] = voronoiDiagram(dt);
    voronoi_con_net = voronoi_contact(cAll, Ncell);
    
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_next_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_next_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX_n = [];
    voronoiY_n = [];
    for i = -1:1
        for j = -1:1
            voronoiX_n = [voronoiX_n,xcomp + lengthscale(end-1) * i];
            voronoiY_n = [voronoiY_n,ycomp + lengthscale(end)* j];
        end
    end

    dt = delaunayTriangulation(voronoiX_n',voronoiY_n');
    [vAll_next,cAll_next] = voronoiDiagram(dt);
    voronoi_con_net_next = voronoi_contact(cAll_next, Ncell);
    
    cell_list = [];
%     for i = (1 + 4 * Ncell) : (5 * Ncell)
%         if length(cAll{i})~=length(cAll_next{i})
%             count = count + 1;
%             cell_list = [cell_list, mod(i-1,Ncell)];
%         end
%     end

    % get cells involved in T1 
    for i = 1 : Ncell
        if ~isequal(voronoi_con_net{i}, voronoi_con_net_next{i})
            count = count + 1;
            cell_list = [cell_list, mod(i-1,Ncell)];
        end
    end
   
    if (count >0)
        % record the contact difference for each cell
        differ_list = {};
        for i = 1 : length(cell_list)
            diff = setdiff(voronoi_con_net{cell_list(i)+1},voronoi_con_net_next{cell_list(i)+1});
            diff = [diff, setdiff(voronoi_con_net_next{cell_list(i)+1},voronoi_con_net{cell_list(i)+1})];
            differ_list{i} =  diff-1; 
        end

        T1_list = [];
        for i = 1 : length(cell_list)
            if length(differ_list{i}) ==1
                T1_sub_list = [cell_list(i)];
                for j = 1 : length(cell_list)
                    % if j is connected to i and the cell which i lost
                    % contact with, then it must be part of T1
                    if ismember(cell_list(j)+1,voronoi_con_net_next{cell_list(i)+1})&&ismember(cell_list(j)+1,voronoi_con_net_next{differ_list{i}+1})
                        T1_sub_list = [T1_sub_list, cell_list(j)];
                    end
                end
                T1_sub_list = [T1_sub_list,differ_list{i}];
                T1_sub_list = unique(T1_sub_list);
                if length(T1_sub_list) == 4
                    T1_list = [T1_list; T1_sub_list];
                end
                
            end
        end

        T1_list = unique(T1_list,'rows');
        count = size(T1_list,1);

        for i = 1: size(T1_list,1)
            if isempty(T1_index)
                T1_cells = [T1_cells; T1_list(i,:)];
                T1_index = [T1_index; frame_i+1];
            else
                % if this is a reverse T1
                [logi,loca]= ismember(T1_list(i,:), T1_cells, 'row');
                if logi
                    count = count - 1;
                    % if recording T1 rate
                    if mode==1
                      T1_count(T1_index(loca)) = T1_count(T1_index(loca)) - 1;
                    end
                    T1_cells(loca,:) = [];
                    T1_index(loca,:) = [];
                else
                    T1_cells = [T1_cells; T1_list(i,:)];
                    T1_index = [T1_index; frame_i+1];
                end
            end
        end
     end
%     else
%      if mod(count,4) > 0
%         disp(cell_list);
%         text_range = (1 + 4 * Ncell) : (5 * Ncell);
%         figure(fig); hold on
%         clf;
%         voronoi(voronoiX,voronoiY)
%         plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell))}, text_range');
%          text(voronoiX(text_range), voronoiY(text_range), plabels, 'FontWeight', ...
%                'bold', 'HorizontalAlignment','center', ...
%                'BackgroundColor', 'none')
%         xlim([0,lengthscale(end-1)]);
%         ylim([0,lengthscale(end)]);
%         hold off
%         
%         figure(fig+1); hold on
%         clf;
%         voronoi(voronoiX_n,voronoiY_n)
%         plabels = arrayfun(@(n) {sprintf('C%d', mod(n-1,Ncell))}, text_range');
%          text(voronoiX_n(text_range), voronoiY_n(text_range), plabels, 'FontWeight', ...
%                'bold', 'HorizontalAlignment','center', ...
%                'BackgroundColor', 'none')
%         xlim([0,lengthscale(end-1)]);
%         ylim([0,lengthscale(end)]);
%         hold off
%      end
%     end
    % if recording number of cumulated T1
    if mode==2
        count = size(T1_cells,1);
    % if recording fraction of cells involved
    elseif mode ==3
        count = size(unique(T1_cells),1)/Ncell;
    end
end


function count = T1_swap_3(xpos_at_frame, ypos_at_frame, xpos_at_next_frame, ypos_at_next_frame, Ncell, lengthscale, frame_i, fig)
    global T1_cells;
    global T1_index;
    global T1_count;
    count = 0;
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX = [];
    voronoiY = [];
    for i = -1:1
        for j = -1:1
            voronoiX = [voronoiX,xcomp + lengthscale(end-1) * i];
            voronoiY = [voronoiY,ycomp + lengthscale(end)* j];
        end
    end

    dt = delaunayTriangulation(voronoiX',voronoiY');
    [vAll,cAll] = voronoiDiagram(dt);
    voronoi_con_net = voronoi_contact(cAll, Ncell);
    
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_next_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_next_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = mod(cx_tmp,lengthscale(end-1));
        ycomp(ci) = mod(cy_tmp,lengthscale(end));
    end
    voronoiX_n = [];
    voronoiY_n = [];
    for i = -1:1
        for j = -1:1
            voronoiX_n = [voronoiX_n,xcomp + lengthscale(end-1) * i];
            voronoiY_n = [voronoiY_n,ycomp + lengthscale(end)* j];
        end
    end

    dt = delaunayTriangulation(voronoiX_n',voronoiY_n');
    [vAll_next,cAll_next] = voronoiDiagram(dt);
    voronoi_con_net_next = voronoi_contact(cAll_next, Ncell);
    
    cell_list = [];
%     for i = (1 + 4 * Ncell) : (5 * Ncell)
%         if length(cAll{i})~=length(cAll_next{i})
%             count = count + 1;
%             cell_list = [cell_list, mod(i-1,Ncell)];
%         end
%     end

    % get cells involved in T1 
    for i = 1 : Ncell
        if ~isequal(voronoi_con_net{i}, voronoi_con_net_next{i})
            count = count + 1;
            cell_list = [cell_list, mod(i-1,Ncell)];
        end
    end
   
    if (count >0)
        % record the contact difference for each cell
        differ_list = {};
        for i = 1 : length(cell_list)
            diff = setdiff(voronoi_con_net{cell_list(i)+1},voronoi_con_net_next{cell_list(i)+1});
            diff = [diff, setdiff(voronoi_con_net_next{cell_list(i)+1},voronoi_con_net{cell_list(i)+1})];
            differ_list{i} =  diff-1; 
        end

        T1_list = [];
        for i = 1 : length(cell_list)
            if length(differ_list{i}) ==1
                T1_sub_list = [cell_list(i)];
                for j = 1 : length(cell_list)
                    % if j is connected to i and the cell which i lost
                    % contact with, then it must be part of T1
                    if ismember(cell_list(j)+1,voronoi_con_net_next{cell_list(i)+1})&&ismember(cell_list(j)+1,voronoi_con_net_next{differ_list{i}+1})
                        T1_sub_list = [T1_sub_list, cell_list(j)];
                    end
                end
                T1_sub_list = [T1_sub_list,differ_list{i}];
                T1_sub_list = unique(T1_sub_list);
                if length(T1_sub_list) == 4
                    T1_list = [T1_list; T1_sub_list];
                end
                
            end
        end

        T1_list = unique(T1_list,'rows');
        count = size(T1_list,1);

        for i = 1: size(T1_list,1)
            if isempty(T1_index)
                T1_cells = [T1_cells; T1_list(i,:)];
                T1_index = [T1_index; frame_i+1];
            else
                % if this is a reverse T1
                [logi,loca]= ismember(T1_list(i,:), T1_cells, 'row');
                if logi
                    count = count - 1;
                    %T1_count(T1_index(loca)) = T1_count(T1_index(loca)) - 1;
                    T1_cells(loca,:) = [];
                    T1_index(loca,:) = [];
                else
                    T1_cells = [T1_cells; T1_list(i,:)];
                    T1_index = [T1_index; frame_i+1];
                end
            end
        end
    end
    count = size(T1_cells,1);
end

function voronoi_con_net = voronoi_contact(cAll, Ncell)
    voronoi_con_net = {};
    for i = (1 + 4 * Ncell) : (5 * Ncell)
        mod_i = mod(i-1,Ncell)+1;
        neighbers = [];
        for element =  cAll{i}
            for j = 1:Ncell*9
                mod_j = mod(j-1,Ncell)+1;
                element_of_other_cell = cAll{j};
                if (sum(element_of_other_cell==element)>0 && mod_j ~= mod_i)
                    neighbers = [neighbers, mod_j];
                end
            end
        end
        voronoi_con_net{mod_i} = sort(unique(neighbers));
    end
end

function order = cal_order(coordinate,N)

    vel = coordinate(1:end - N,:) - coordinate(N+1:end,:);
    velx = reshape(vel(:,1) , N, []);
    vely = reshape(vel(:,2) , N, []);
    vel_norm = sqrt(velx.^2+ vely.^2);
    velx = mean(velx./vel_norm,1);
    vely = mean(vely./vel_norm,1);
    
%     velx = mean(velx,1);
%     vely = mean(vely,1);
    
    order = sqrt(velx.^2+ vely.^2);
    order = mean(order,'all');
      
end

function order = cal_order1(vel)

    order = sqrt(vel(:,1).^2+ vel(:,2).^2);
    order = mean(order,'all');
      
end




function isequal = compare_contact(voronoi_con_net,jammed_voronoi_con_net,Ncell)
    isequal = 1;
    for i = 1:Ncell
        if sum(ismember(jammed_voronoi_con_net{i},voronoi_con_net{i}),'all') < 3
            isequal = 0;
        end
    end
end

function isequal = compare_contact1(voronoi_con_net,jammed_voronoi_con_net,xcomp,ycomp,lengthscale,Ncell)
    isequal = 1;
    for i = 1:Ncell
        changedcontact = find(~ismember(voronoi_con_net{i},jammed_voronoi_con_net{i}));
        for index = changedcontact
            dis1 = cell_distance(i,voronoi_con_net{i}(index),xcomp,ycomp,lengthscale);
            dis2 = mean(cell_distance(i,jammed_voronoi_con_net{i},xcomp,ycomp,lengthscale),'all');
            if dis1> 2 * dis2
                isequal = 0;
            end
        end
    end
end

function distance = cell_distance(i,j,xcomp,ycomp,lengthscale)
    ci_x = xcomp(i);
    ci_y = ycomp(i);
    cj_x = xcomp(j);
    cj_y = ycomp(j);
    dx = ci_x - cj_x;
    dy = ci_y - cj_y;
    dx = dx - lengthscale(end-1)*round(dx./lengthscale(end-1));
    dy = dy - lengthscale(end)*round(dy./lengthscale(end));
    distance = sqrt(dx.^2+dy.^2);
end

function n_particle = density_fluctuation(N, coordinate, frames, Ncell, lengthscale)
    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);
    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
    end
    
    xcomp = mod(xcomp,lengthscale(end-1));
    ycomp = mod(ycomp,lengthscale(end));
    n_particle = xcomp<3*lengthscale(end-1)/6 &  xcomp>lengthscale(end-1)/6 & ycomp<3*lengthscale(end)/6 & ycomp>lengthscale(end)/6;
    n_particle = sum(n_particle,2);
    
end



function [xcomp,ycomp]=cal_c_pos(xpos_at_frame, ypos_at_frame, Ncell, lengthscale)
    xcomp = zeros(1,Ncell);
    ycomp = zeros(1,Ncell);
    for ci = 1:Ncell
        index = sum(lengthscale(1:ci),'all');
        start_point_last = 1 + index - lengthscale(ci);
        end_point_last = index;
        cx_tmp = mean(xpos_at_frame(start_point_last:end_point_last),'all');
        cy_tmp = mean(ypos_at_frame(start_point_last:end_point_last),'all');
        xcomp(ci) = cx_tmp;
        ycomp(ci) = cy_tmp;
    end
end


function [MSD,deltaT] = cal_msd(N, coordinate, frames, Ncell, lengthscale)


    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);
    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
    end
    
    
    % create MSD array (y-axis of MSD plot)
    MSD = zeros(round(9*frames/10),1);
    NT = length(MSD);
    % loop over the different possible time windows, calculate MSD for each
    % time window size
    for ii = 1:NT

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

function [ISF,deltaT] = cal_ISF(N, coordinate, frames, Ncell, lengthscale)


    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);
    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
    end
    
    
    % create MSD array (y-axis of MSD plot)
    ISF = zeros(round(9*frames/10),1);
    NT = length(ISF);
    q = sqrt(lengthscale(end)*lengthscale(end-1)/(3.14*Ncell));
    logindex = unique(round(logspace(0, log10(NT),100)));
    % loop over the different possible time windows, calculate MSD for each
    % time window size
    for ii = logindex

        % calculate x displacements, separated by ii indices
        dx = xcomp(1+ii:end,:) - xcomp(1:end-ii,:);

        % calculate y displacements similarly
        dy = ycomp(1+ii:end,:) - ycomp(1:end-ii,:);
    
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

function [structure, wave] = cal_structure_factor(N, coordinate, frames, Ncell, lengthscale)


    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);
    for i = 1 :1
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
    end
    
    
    % create MSD array (y-axis of MSD plot)
    structure = [];
    wave = [];
    r0 = sqrt(0.92 * lengthscale(end)*lengthscale(end-1)/(3.14*Ncell));

    for q = 2*pi/lengthscale(end) : 0.1 : 2*pi/(0.5*r0)

        % calculate x displacements, separated by ii indices
        dx = xcomp(1,:);

        % calculate y displacements similarly
        dy = ycomp(1,:);
    
        count = 0;
        s = 0;
        for th = 0: 0.1: 2*3.141
        % take mean over all displacements
            s = s + abs(sum(exp(sqrt(-1) * (cos(th)*q*dx + sin(th)*q*dy)),'all'))^2;
            count = count + 1;
        end
        % store in MSD array
        structure = [structure, s/(count* Ncell)];
        wave = [wave, q*r0/pi];
    end


end

function [std_n, mean_n] = giant_number(N, coordinate, frames, Ncell, lengthscale)
    xcomp=zeros(frames,Ncell);
    ycomp=zeros(frames,Ncell);

    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        [xcomp_t,ycomp_t]=cal_c_pos(coordinate(start_point:end_point,1),coordinate(start_point:end_point,2), Ncell, lengthscale);
        xcomp(i,:)= xcomp_t;
        ycomp(i,:)= ycomp_t;
 
    end
    
    xcomp = mod(xcomp,lengthscale(end-1));
    ycomp = mod(ycomp,lengthscale(end));
    
    mean_n = [];
    std_n = [];
    for i = 1:5
        %n_particle = xcomp<i*lengthscale(end-1)/10 &  xcomp>0 & ycomp<i*lengthscale(end)/10 & ycomp>0;
        n_particle = xcomp<i*lengthscale(end-1)/10 &  xcomp>0 & ycomp<(1/2 + i/20)*lengthscale(end) & ycomp>(1/2 - i/20)*lengthscale(end);
        n_particle = sum(n_particle,2);
        std_n = [std_n, std(n_particle)];
        mean_n = [mean_n, mean(n_particle)];
    end
end


function [breaked, new] = compare_contact2(voronoi_con_net,jammed_voronoi_con_net,Ncell)
    breaked = 0;
    new = 0;
    for i = 1:Ncell
        breaked = breaked + length(jammed_voronoi_con_net{i}) - sum(ismember(voronoi_con_net{i}, jammed_voronoi_con_net{i}),'all');
        new = new + length(voronoi_con_net{i}) - sum(ismember(voronoi_con_net{i}, jammed_voronoi_con_net{i}),'all');
    end
end



function [MSD,deltaT] = cal_msd_vertex(N, coordinate, frames)


    xcomp=zeros(frames,N);
    ycomp=zeros(frames,N);
    for i = 1 :frames
        start_point = 1 + N * ( i - 1 );
        end_point = N * i;
        xcomp(i,:)= coordinate(start_point:end_point,1);
        ycomp(i,:)= coordinate(start_point:end_point,2);
    end
    
    
    % create MSD array (y-axis of MSD plot)
    MSD = zeros(round(9*frames/10),1);
    NT = length(MSD);
    % loop over the different possible time windows, calculate MSD for each
    % time window size
    for ii = 1:NT

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






% kb=sort(v0(:,3));
% figure(7)
% plot(kb,all_mean_cal_A)
% xlabel('kb');ylabel('calA');

% for ground shape
% kb = unique(v0(:,3));
% kl = unique(v0(:,4));
% P = polyfit(log10(kl'),log10(1.12-all_mean_cal_A(3,:)),1)
% loglog(kl,1.12-all_mean_cal_A(3,:))
% P = polyfit(log10(kb(2:end)),log10(1.12-all_mean_cal_A(2:end,3)),1)
% loglog(kb,1.12-all_mean_cal_A(:,3))

% figure(7);hold on
% offset = 1.08;
% scatter(all_mean_cal_A(:)-offset,ifjammed(:))
% P = polyfit(log10(all_mean_cal_A(:)-offset), log10(ifjammed(:)), 1);
% yfit = P(1)*log10(all_mean_cal_A(:)-offset)+P(2);
% plot(all_mean_cal_A(:)-offset,10.^(yfit),'r-.');
% ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
% xlabel('calA_m-calA_g');ylabel('MSD');

% figure(7);hold on
% scatter(all_mean_cal_A(:),ifjammed(:))
% P = polyfit(log10(all_mean_cal_A(:)), log10(ifjammed(:)), 1);
% yfit = P(1)*log10(all_mean_cal_A(:))+P(2);
% plot(all_mean_cal_A(:),10.^(yfit),'r-.');
% ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
% xlabel('calA');ylabel('MSD');


