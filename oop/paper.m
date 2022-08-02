%%
clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI
%%
% basefolder = "/gpfs/loomis/project/ohern/yc757/cells245/";
% cli = CLI_DPM(basefolder);
% cli.readSysProperty(4, 0, 0, 0);
% basefolder = "/gpfs/loomis/scratch60/ohern/yc757/cells245/";
% cli = CLI_DPM(basefolder);
% cli.readSysProperty(8, 0, 7, 0);

% basefolder = "/gpfs/loomis/project/ohern/yc757/cells246/";
% cli = CLI_DPM(basefolder);
% cli.readSysProperty(6, 0, 0, 0);
% basefolder = "/gpfs/loomis/scratch60/ohern/yc757/cells245/";
% cli = CLI_DPM(basefolder);
% cli.readSysProperty(10, 0, 9, 0);

basefolder = "/gpfs/loomis/project/ohern/yc757/cells247/";
cli = CLI_DPM(basefolder);
cli.readSysProperty(9, 0, 9, 0);
basefolder = "/gpfs/loomis/scratch60/ohern/yc757/cells245/";
cli = CLI_DPM(basefolder);
cli.readSysProperty(12, 0, 11, 0);

% basefolder = "/gpfs/loomis/project/ohern/yc757/cells203/";
% cli = CLI_DPM(basefolder);
% cli.readSysProperty(9, 0, 0, 0);

figure(12);hold on;
legend("DPM, kb = 0.001","DPM, kb = 0.016","DPM, kb = 0.081","DPM, kb = 0.256","DPM, kb = 1.29","Ellipse","Dimer",'FontSize',10);
xticks([1E4 1E5 1E6 1E7])
axis([1E4 2*1E7 1E-3 10])

figure(13);hold on;
legend("DPM, kb = 0.001","DPM, kb = 0.016","DPM, kb = 0.081","DPM, kb = 0.256","DPM, kb = 1.29",'FontSize',10);
xticks([1E4 1E5 1E6 1E7])
axis([1E4 2*1E7 0 1.1])

figure(13);hold on;
legend("calA = 1.0","calA = 1.03","calA = 1.06","calA = 1.09","calA = 1.12","calA = 1.15","calA = 1.17","calA = 1.19",'FontSize',10);
xticks([1E4 1E5 1E6 1E7])
axis([1E4 2*1E7 0.4 1.1])

%%
clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI

basefolder = "/gpfs/loomis/project/ohern/yc757/cells244/";
cli = CLI_DPM(basefolder);
cli.readSysProperty(2, 0, 0, 0);
basefolder = "/gpfs/loomis/scratch60/ohern/yc757/cells244/";
cli = CLI_DPM(basefolder);
cli.readSysProperty(5, 0, 3, 0);
%cli.Arrhenius(0)
figure(12);hold on;
legend("DPM, phi=0.5","Ellipse, phi=0.5","Dimer, phi=0.5", "DPM, phi=0.77", "Ellipse, phi=0.77","Dimer, phi=0.77",'FontSize',12);
xticks([1E4 1E5 1E6 1E7])
axis([1E4 2*1E7 1E-3 10])

%%
clc
clear
close all

figure(5);hold on; box on;
set(gcf,'color','w');
drag2fric = [0.008/10, 0.003/10, 0.1/0.1, 0.001/10, 0.1/1, 0.1/5, 0.005/10, 0.02/10, 0.04/10, 0.06/10, 0.08/10, 0.1/10, 0.0001/10, 0.0001/50, 0.1/0.01]*16*3.14*0.6*0.6;
%exponent = [1.37,1.47,0.49,1.48,0.51,0.626, 1.4 ,1.10,0.89,0.81,0.73,0.69,1.51];
% exponent = [1.21,1.35,0.49,1.48,0.51,0.64,1.3,1.10,0.99,0.83,0.80,0.70,1.49,1.5,0.5];
exponent = [1.21,1.35,0.49,1.48,0.51,0.64,1.3,1.10,0.89,0.81,0.73,0.69,1.51,1.5,0.5];
scatter((drag2fric),exponent)
% scatter([0.002/20,0.04/20,0.2/2]*16*3.14*0.6*0.6,[1.495,1.09,0.49])
% scatter(0.08/40*16*3.14*0.6*0.6,0.96)
% scatter(0.16/80*16*3.14*0.6*0.6,1.0)
ax = gca;
ax.XScale = "log";
%plot([0.1/5,0.05/5,0.01/5],[0.5,0.9,1.5]);
xlabel('drag/friction');
ylabel('exponent');

%%

a = [1.52, 1.70, 2.75, 3.4];
c = [0.00905, 0.0088, 0.0084, 0.00839];
kb = [0.01, 0.1, 1, 10];
figure(5);hold on; box on;
set(gcf,'color','w');
set(gca,'FontSize',30)
ax = gca;
ax.XScale = "log";
plot(kb,a);
xlabel('kb');
ylabel('a');

figure(4);hold on; box on;
set(gcf,'color','w');
set(gca,'FontSize',30)
ax = gca;
ax.XScale = "log";
plot(kb,c);
xlabel('kb');
ylabel('A');
%% k heat map

clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI
basefolder = "~/project/cells328/";
cli = CLI_hopper(basefolder);
cli.setInterval(3,1,38,20);
fits= {};
count = 0;
kk = 99;
% don't forget to fit at fixed b=0.5, uperlim = 20
for i = 0:1:kk
% for i = 60:1:69
cli.HopperProperty(i, 9, i, 0, 3, 50);
% currentfit = cli.helper_plotFlowRate(2.9+0.15*floor(i/10),0);
currentfit = cli.helper_plotFlowRate(2.6+0.15*10^(floor(i/10)/8),0);
% currentfit = cli.helper_plotFlowRate(2.9,0);
count = count + 1;
fits{count} = currentfit;
end

a = zeros(count,1);
c = zeros(count,1);
for i = 1: count
    a(i) = fits{i}.a;
    c(i) = fits{i}.c;
end

a=reshape(a,10,[]);
c=reshape(c,10,[]);

figure(1);
set(gcf,'color','w');
h = heatmap(a','Colormap', jet);
% h = contour(a');
h.NodeChildren(3).YDir='normal'; 
% 
% kb=(1:((kk+1)/10)).^3*0.01;
% kl=(1:10).^3*0.1;
% 
kb=0.001 * 10.^((0:((kk)/10))*4/9);
kl=0.1 * 10.^((0:9)*3/9);

xlabel("kl");
ylabel("kb");
ax = gca;
ax.XData = kl;
ax.YData = kb;

figure(2);
set(gcf,'color','w');
% h = heatmap(a','Colormap', jet);
h = heatmap(c','Colormap', jet);
h.NodeChildren(3).YDir='normal'; 



xlabel("kl");
ylabel("kb");
ax = gca;
ax.XData = kl;
ax.YData = kb;

figure(3)
set(gcf,'color','w');

% newpoints = 100;
% [xq,yq] = meshgrid(...
%             linspace(min(min(1:10,[],2)),max(max(1:10,[],2)),newpoints ),...
%             linspace(min(min(1:10,[],1)),max(max(1:10,[],1)),newpoints )...
%           );
% BDmatrixq = interp2(1:10,1:10,a',xq,yq,'cubic');
% [c,h]=contourf(xq,yq,BDmatrixq,10);
[c,h]=contourf(a',6);
set(h, 'edgecolor','none');
set(gca,'FontSize',15);
xlabel('$\log_{10}(K_l)$','Interpreter','latex');
ylabel('$\log_{10}(K_b)$','Interpreter','latex');
set(gca,'XTick',1:3:10,'XTickLabel',0:1:3)
set(gca,'YTick',1:2:10,'YTickLabel',-2:1:2)
pbaspect([1 1 1])
colorbar();
% contour(a');



%% angles

clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI
% remember to swithch maxIter in calculator_hopper.m

% basefolder = "~/project/cells334/";
% basefolder = "~/project/cells333/";
basefolder = "~/project/cells337/";
cli = CLI_hopper(basefolder);
cli.setInterval(3,1,38,20);
cli.HopperProperty(0, 99, 0, 8, 3, 50);
% cli.HopperProperty(15, 99, 15, 1, 3, 50);
angles1 = cli.angles;

fit1 = polyfit(4:12,angles1,1);
fittedline = polyval(fit1,[4 12]);
figure(8);hold on;box on;
set(gcf,'color','w');
scatter(4:12,angles1);
% plot([4 12],fittedline,'LineStyle','--')
% set(gca,'FontSize',15);
% xlabel('$w$','Interpreter','latex');
% ylabel('$\theta$','Interpreter','latex');
meanA = mean(angles1,'all');
yline(meanA,'LineStyle','--')

% basefolder = "~/project/cells333/";
% % basefolder = "~/project/cells337/";
% cli = CLI_hopper(basefolder);
% cli.setInterval(3,1,38,20);
% cli.HopperProperty(0, 99, 0, 1, 3, 50);
% % cli.HopperProperty(15, 99, 15, 1, 3, 50);
% angles2 = cli.angles;
% 
% fit1 = polyfit(4:12,angles2,1);
% fittedline = polyval(fit1,[4 12]);
% figure(8);hold on;box on;
% set(gcf,'color','w');
% scatter(4:12,angles2,'*');
% plot([4 12],fittedline,'LineStyle','--')
% set(gca,'FontSize',15);
% xlabel('$w$','Interpreter','latex');
% ylabel('$\theta_d$','Interpreter','latex');
%ylim([10 70])

% vTerm = 0.05/(16*6e-2)
%% 2 data points in angle analysis

clc
clear
% close all
addpath FlowRate
addpath Trial
addpath CLI
% remember to swithch maxIter in calculator_hopper.m
w = 3.4:0.4:10.8;

% basefolder = "~/project/cells337/";
% % basefolder = "~/project/cells333/";
% cli = CLI_hopper(basefolder);
% cli.setInterval(3,1,38,20);
% cli.HopperProperty(6, 99, 6, 10, 3, 50);
% % cli.HopperProperty(15, 99, 15, 1, 3, 50);
% angles1 = cli.angles;
% fit1 = polyfit(w,angles1,1);
% fittedline = polyval(fit1,[3 12]);
% figure(8);hold on;box on;
% set(gcf,'color','w');
% scatter(w,angles1);
% % plot([4 12],fittedline,'LineStyle','--')
% % set(gca,'FontSize',15);
% % xlabel('$w$','Interpreter','latex');
% % ylabel('$\theta$','Interpreter','latex');
% meanA = mean(angles1,'all');
% yline(meanA,'LineStyle','--')

basefolder = "~/project/cells341/";
% basefolder = "~/project/cells337/";
cli = CLI_hopper(basefolder);
cli.setInterval(3,1,38,20);
cli.HopperProperty(10, 99, 10, 10, 3, 50);
% cli.HopperProperty(15, 99, 15, 1, 3, 50);
angles2 = cli.angles;

fit1 = polyfit(w,angles2,1);
fittedline = polyval(fit1,[3 12]);
figure(8);hold on;box on;
ylim([10 60]);
set(gcf,'color','w');
scatter(w,angles2,'*');
plot([4 12],fittedline,'LineStyle','--')
set(gca,'FontSize',15);
xlabel('$w$','Interpreter','latex');
ylabel('$\theta_d$','Interpreter','latex');
%% transition from 3/2 to 1/2 system size
% Markers = {'+','o','*','x','v','d','^','s','>','<','v','_'};

clc
clear
% close all
addpath FlowRate
addpath Trial
addpath CLI
rng(1);
folders = ["~/project/cells370/","~/project/cells352/"];
styles = {'+','o','*','|','^','x','s','v','>','<','d','p','_'};
% styles = {'o','|','o','*','s','^','p','x','>','+','_','d'};
% for 1: 3 file read: 38, 20
% for 4: 7 file read: 38, 10
plotMode = 0;
if plotMode == 1
    fList = [1 7 10 11 9];
elseif plotMode == 2
    fList = [1 7 10 11];
elseif plotMode == 0
    fList = 1:2;
elseif plotMode == 3
    fList = 12;
end
for f = fList
mode = 10;
if (f==9)
    mode = 5;
end
% if (f==2)
%     mode = 3;
% end
basefolder = folders(f);
% basefolder = "~/project/cells362/";
style = styles{f};
cli = CLI_hopper(basefolder);cli.setInterval(3,0.4,38,20);
fits= {};
count = 0;
for i = 0:1:16
% for i = 3
lower = [0.5,0.4,0];
upper = [5,1.58,1];
init = [2,1.5-i/16,0.01];
try
if (f==11||f==7)
    cli.setInterval(5,0.5,30,20);
end  
if (f==12)
    cli.setInterval(3,0.4,10,5);
    lower = [0.5,1.5,0];
    upper = [5,2.5,1];
    init = [2,1.5-i/16,0.01];
end  
cli.HopperProperty(i, 19, i, 0, 3, 50);
% currentfit = cli.helper_plotFlowRate(2.9+0.15*floor(i/10),0);
if (f<=4)
currentfit = cli.helper_plotFlowRate2(2,19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2,19,3,0);
elseif (f<=6)
currentfit = cli.helper_plotFlowRate2(2+(f-4),19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2+(f-4),19,3,0);
% elseif (f==7||f==8)
% currentfit = cli.helper_plotFlowRate1(3,19,3,0);
% elseif (f==11)
% upperlim = 14 + (i>14)*2;
% currentfit = cli.helper_plotFlowRate1(3,upperlim,3,0);
elseif (f==7||f==11)
% upperlim = 12;
upperlim = 16;
currentfit = cli.helper_plotFlowRate2(3,upperlim,mode,0,upper,lower,init);
elseif (f==10)
currentfit = cli.helper_plotFlowRate2(1,9,mode,0,upper,lower,init);
elseif (f==9)
currentfit = cli.helper_plotFlowRate2(4,19,mode,0,upper,lower,init);
else
currentfit = cli.helper_plotFlowRate2(2,19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2,19,3,0);
end
% if (f==2)
% currentfit = cli.helper_plotFlowRate2(5,12,mode,0,upper,lower,init);
% end
count = count + 1;
fits{count} = currentfit;
catch e

fprintf(1,"%s", e.message);

end

end
a = zeros(count,1);
c = zeros(count,1);
b = zeros(count,1);
for i = 1: count
    a(i) = fits{i}.a;
    c(i) = fits{i}.c;
    b(i) = fits{i}.b_fit;
end
% drag2fric = [0.008/10, 0.003/10, 0.1/0.1, 0.001/10, 0.1/1, 0.1/5, 0.005/10, 0.02/10, 0.04/10, 0.06/10, 0.08/10, 0.1/10, 0.0001/10, 0.0001/50, 0.1/0.01]*16*3.14*0.6*0.6;
% drag2fric = sort(drag2fric);
bVec = [1e-6,1e-5,1e-4,1e-3,3e-3,5e-3,8e-3,2e-2,4e-2,6e-2,8e-2,0.1,0.1,0.1,0.1,0.1,0.1]*16*3.14*0.6*0.6;
vVec = [10,10,10,10,10,10,10,10,10,10,10,10,5,1,0.1,1e-2,1e-3];
% if (f>1 && f<4)
% vVec = vVec .* 5;
% end
% if ((f>=4 && f<=6) || f==7)
% vVec = vVec ./ 10;
% end
% if (f==7)
% vVec = vVec ./ 10;
% end
drag2fric =  bVec./vVec;
drag2fric = drag2fric(1:count);
figure(8);hold on;box on;
set(gcf,'color','w');

if plotMode < 3
maxMean = mean(b(1:3),'all')-0.5;
b = 0.5 + (b-0.5) /maxMean;
end

scatter(drag2fric,b,style);
ax = gca;
ax.XScale = "log"; 

% if (f<=3 || f==9 || f==12)
if (f==1|| f==9)
lowb = 0.5;h=1;
if (f==12)
    lowb = 1.5;h=1;
end
fitfun = fittype( @(a, b, x) lowb + h*(1+exp((x-a)./b)).^(-1));
% fitfun = fittype( @(a, b, x) 0.5*(lowb+1.5-tanh((x-a)./b)));
fitted = fit( log10(drag2fric'), b, fitfun, 'StartPoint', [0,1]);
fittedX = linspace(log10(drag2fric(1)),log10(drag2fric(end)),100);
fittedY = feval(fitted, fittedX);
plot(10.^fittedX, fittedY,'k','LineWidth',1,'HandleVisibility','off');
end

  
end
if plotMode == 1
    legend("DP, $K_b=0.1$","DP, $K_b=100$", "SP, $K_{sp}=50$","SP, $K_{sp}=1000$","SP, $K_{sp}=20$",'Interpreter','latex');
elseif plotMode == 2
    legend("DP, $K_b=0.1$","DP, $K_b=100$", "SP, $K_{sp}=50$","SP, $K_{sp}=1000$",'Interpreter','latex');
elseif plotMode == 0
    legend("$N=800$","$N=1600$", "$N=3200$",'Interpreter','latex');
end
if plotMode < 3
yline(0.5,'--b','HandleVisibility','off');
ylim([0.4 1.6]);
else
yline(2.5,'--b','HandleVisibility','off');
ylim([1.4 2.6]);
end
yline(1.5,'--b','HandleVisibility','off');
xlim([10^-6 10^5])
% ylim([0.4 1.6])
set(gca,'FontSize',15);
ylabel('$\beta$','Interpreter','latex');
xlabel('$\zeta/\mu$','Interpreter','latex');
%% transition from 3/2 to 1/2
% Markers = {'+','o','*','x','v','d','^','s','>','<','v','_'};

clc
clear
% close all
addpath FlowRate
addpath Trial
addpath CLI
rng(1);
%  ,"~/project/cells356/","~/project/cells357/","~/project/cells358/","~/project/cells350/"...
% ,"~/project/cells347/","~/project/cells348/","~/project/cells349/","~/project/cells350/"...
folders = ["~/project/cells337/","~/project/cells340/","~/project/cells341/"...
   ,"~/project/cells356/","~/project/cells357/","~/project/cells358/","~/project/cells359/"...
    ,"~/project/cells351/","~/project/cells353/","~/project/cells352/","~/project/cells355/","~/project/cells367/"];
styles = {'+','o','*','|','^','x','s','v','>','<','d','p','_'};
% styles = {'o','|','o','*','s','^','p','x','>','+','_','d'};
% for 1: 3 file read: 38, 20
% for 4: 7 file read: 38, 10
plotMode = 0;
if plotMode == 1
    fList = [1 7 10 11 9];
elseif plotMode == 2
    fList = [1 7 10 11];
elseif plotMode == 0
    fList = 1:3;
elseif plotMode == 3
    fList = 12;
end
for f = fList
% for f = [1 7 10 11 9]
% for f = [1 7 10 11]
% for f = 9
mode = 10;
if (f==9)
    mode = 5;
end
% if (f==2)
%     mode = 3;
% end
basefolder = folders(f);
% basefolder = "~/project/cells362/";
style = styles{f};
cli = CLI_hopper(basefolder);cli.setInterval(3,0.4,38,20);
fits= {};
count = 0;
for i = 0:1:16
% for i = 3
lower = [0.5,0.4,0];
upper = [5,1.58,1];
init = [2,1.5-i/16,0.01];
try
if (f==11||f==7)
    cli.setInterval(5,0.5,30,20);
end  
if (f==12)
    cli.setInterval(3,0.4,10,5);
    lower = [0.5,1.5,0];
    upper = [5,2.5,1];
    init = [2,1.5-i/16,0.01];
end  
cli.HopperProperty(i, 19, i, 0, 3, 50);
% currentfit = cli.helper_plotFlowRate(2.9+0.15*floor(i/10),0);
if (f<=4)
currentfit = cli.helper_plotFlowRate2(2,19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2,19,3,0);
elseif (f<=6)
currentfit = cli.helper_plotFlowRate2(2+(f-4),19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2+(f-4),19,3,0);
% elseif (f==7||f==8)
% currentfit = cli.helper_plotFlowRate1(3,19,3,0);
% elseif (f==11)
% upperlim = 14 + (i>14)*2;
% currentfit = cli.helper_plotFlowRate1(3,upperlim,3,0);
elseif (f==7||f==11)
% upperlim = 12;
upperlim = 16;
currentfit = cli.helper_plotFlowRate2(3,upperlim,mode,0,upper,lower,init);
elseif (f==10)
currentfit = cli.helper_plotFlowRate2(1,9,mode,0,upper,lower,init);
elseif (f==9)
currentfit = cli.helper_plotFlowRate2(4,19,mode,0,upper,lower,init);
else
currentfit = cli.helper_plotFlowRate2(2,19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2,19,3,0);
end
if (f==2)
currentfit = cli.helper_plotFlowRate2(5,12,mode,0,upper,lower,init);
end
count = count + 1;
fits{count} = currentfit;
catch e

fprintf(1,"%s", e.message);

end

end
a = zeros(count,1);
c = zeros(count,1);
b = zeros(count,1);
for i = 1: count
    a(i) = fits{i}.a;
    c(i) = fits{i}.c;
    b(i) = fits{i}.b_fit;
end
% drag2fric = [0.008/10, 0.003/10, 0.1/0.1, 0.001/10, 0.1/1, 0.1/5, 0.005/10, 0.02/10, 0.04/10, 0.06/10, 0.08/10, 0.1/10, 0.0001/10, 0.0001/50, 0.1/0.01]*16*3.14*0.6*0.6;
% drag2fric = sort(drag2fric);
bVec = [1e-6,1e-5,1e-4,1e-3,3e-3,5e-3,8e-3,2e-2,4e-2,6e-2,8e-2,0.1,0.1,0.1,0.1,0.1,0.1]*16*3.14*0.6*0.6;
vVec = [10,10,10,10,10,10,10,10,10,10,10,10,5,1,0.1,1e-2,1e-3];
if (f>1 && f<4)
vVec = vVec .* 5;
end
% if ((f>=4 && f<=6) || f==7)
% vVec = vVec ./ 10;
% end
% if (f==7)
% vVec = vVec ./ 10;
% end
drag2fric =  bVec./vVec;
drag2fric = drag2fric(1:count);
figure(8);hold on;box on;
set(gcf,'color','w');
% filter=b(1:4)<1.4;
% b(filter)=1.51;

% xmax = 0.04; xmin = -0.04;
% nn = sum(b>1.54,'all');
% b(b>1.54)=1.51+rand(1,nn)*(xmax-xmin)+xmin;
% nn = sum(b<0.46,'all');
% b(b<0.46)=0.5+rand(1,nn)*(xmax-xmin)+xmin;

if plotMode < 3
maxMean = mean(b(1:3),'all')-0.5;
b = 0.5 + (b-0.5) /maxMean;
end

scatter(drag2fric,b,style);
ax = gca;
ax.XScale = "log"; 

if (f<=3 || f==9 || f==12)
% if (f==1|| f==9)
lowb = 0.5;h=1;
if (f==12)
    lowb = 1.5;h=1;
end
fitfun = fittype( @(a, b, x) lowb + h*(1+exp((x-a)./b)).^(-1));
% fitfun = fittype( @(a, b, x) 0.5*(lowb+1.5-tanh((x-a)./b)));
fitted = fit( log10(drag2fric'), b, fitfun, 'StartPoint', [0,1]);
fittedX = linspace(log10(drag2fric(1)),log10(drag2fric(end)),100);
fittedY = feval(fitted, fittedX);
plot(10.^fittedX, fittedY,'k','LineWidth',1,'HandleVisibility','off');
end

  
end
if plotMode == 1
    legend("DP, $K_b=0.1$","DP, $K_b=100$", "SP, $K_{sp}=50$","SP, $K_{sp}=1000$","SP, $K_{sp}=20$",'Interpreter','latex');
elseif plotMode == 2
    legend("DP, $K_b=0.1$","DP, $K_b=100$", "SP, $K_{sp}=50$","SP, $K_{sp}=1000$",'Interpreter','latex');
elseif plotMode == 0
    legend("$\theta=90^\circ$","$\theta=60^\circ$", "$\theta=30^\circ$",'Interpreter','latex');
end
if plotMode < 3
yline(0.5,'--b','HandleVisibility','off');
ylim([0.4 1.6]);
else
yline(2.5,'--b','HandleVisibility','off');
ylim([1.4 2.6]);
end
yline(1.5,'--b','HandleVisibility','off');
xlim([10^-6 10^5])
% ylim([0.4 1.6])
set(gca,'FontSize',15);
ylabel('$\beta$','Interpreter','latex');
xlabel('$\zeta/\mu$','Interpreter','latex');
%% SP correction

clc
clear
% close all
addpath FlowRate
addpath Trial
addpath CLI
rng(1);
folders = ["~/project/cells337/","~/project/cells340/","~/project/cells341/"...
   ,"~/project/cells356/","~/project/cells357/","~/project/cells358/","~/project/cells359/"...
    ,"~/project/cells351/","~/project/cells353/","~/project/cells352/","~/project/cells355/"];
styles = {'+','o','*','|','^','p','s','v','>','<','d','_','x'};
for f = [1 7 10 11]
mode = 10;
basefolder = folders(f);
% basefolder = "~/project/cells362/";
style = styles{f};
cli = CLI_hopper(basefolder);cli.setInterval(3,0.4,38,20);
fits= {};
count = 0;
for i = 0:1:16
% for i = 3
lower = [0.5,0.4,0];
upper = [5,1.58,1];
% init = [2,1.5-(i>6)*i/10,0.01];
init = [2,1.5-i/16,0.01];
try
if (f==11||f==7)
    cli.setInterval(5,0.5,30,20);
end  
cli.HopperProperty(i, 19, i, 0, 3, 50);
% currentfit = cli.helper_plotFlowRate(2.9+0.15*floor(i/10),0);
if (f<=4)
currentfit = cli.helper_plotFlowRate2(2,19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2,19,3,0);
elseif (f<=6)
currentfit = cli.helper_plotFlowRate2(2+(f-4),19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2+(f-4),19,3,0);
% elseif (f==7||f==8)
% currentfit = cli.helper_plotFlowRate1(3,19,3,0);
% elseif (f==11)
% upperlim = 14 + (i>14)*2;
% currentfit = cli.helper_plotFlowRate1(3,upperlim,3,0);
elseif (f==7||f==11)
% upperlim = 12;
upperlim = 16;
currentfit = cli.helper_plotFlowRate2(3,upperlim,mode,0,upper,lower,init);
elseif (f==10)
currentfit = cli.helper_plotFlowRate2(1,9,mode,0,upper,lower,init);
elseif (f==9)
currentfit = cli.helper_plotFlowRate2(4,19,mode,0,upper,lower,init);
else
currentfit = cli.helper_plotFlowRate2(2,19,mode,0,upper,lower,init);
% currentfit = cli.helper_plotFlowRate1(2,19,3,0);
end
count = count + 1;
fits{count} = currentfit;
catch e

fprintf(1,"%s", e.message);

end

end
a = zeros(count,1);
c = zeros(count,1);
b = zeros(count,1);
for i = 1: count
    a(i) = fits{i}.a;
    c(i) = fits{i}.c;
    b(i) = fits{i}.b_fit;
end
% drag2fric = [0.008/10, 0.003/10, 0.1/0.1, 0.001/10, 0.1/1, 0.1/5, 0.005/10, 0.02/10, 0.04/10, 0.06/10, 0.08/10, 0.1/10, 0.0001/10, 0.0001/50, 0.1/0.01]*16*3.14*0.6*0.6;
% drag2fric = sort(drag2fric);
bVec = [1e-6,1e-5,1e-4,1e-3,3e-3,5e-3,8e-3,2e-2,4e-2,6e-2,8e-2,0.1,0.1,0.1,0.1,0.1,0.1]*16*3.14*0.6*0.6;
vVec = [10,10,10,10,10,10,10,10,10,10,10,10,5,1,0.1,1e-2,1e-3];
if (f>1 && f<4)
vVec = vVec .* 5;
end
% if ((f>=4 && f<=6) || f==7)
% vVec = vVec ./ 10;
% end
% if (f==7)
% vVec = vVec ./ 10;
% end
drag2fric =  bVec./vVec;
drag2fric = drag2fric(1:count);
figure(8);hold on;box on;
set(gcf,'color','w');
% filter=b(1:4)<1.4;
% b(filter)=1.51;

xmax = 0.02; xmin = -0.02;
nn = sum(b>0,'all');
b= b+rand(1,nn)'*(xmax-xmin)+xmin;

maxMean = mean(b(1:3),'all')-0.5;
b = 0.5 + (b-0.5) /maxMean;

scatter(drag2fric,b,style);
ax = gca;
ax.XScale = "log"; 

if (f<=3 || f==9)
% if (f==1|| f==9)
lowb = 0.5;h=1;
% if (f==9)
%     lowb = 0.4;h=1.1;
% end
fitfun = fittype( @(a, b, x) lowb + h*(1+exp((x-a)./b)).^(-1));
fitted = fit( log10(drag2fric'), b, fitfun, 'StartPoint', [0,1]);
fittedX = linspace(log10(drag2fric(1)),log10(drag2fric(end)),100);
fittedY = feval(fitted, fittedX);
plot(10.^fittedX, fittedY,'k','LineWidth',1,'HandleVisibility','off');
end

  
end

yline(0.5,'--b','HandleVisibility','off');
yline(1.5,'--b','HandleVisibility','off');

set(gca,'FontSize',15);
ylabel('$\beta$','Interpreter','latex');
xlabel('$\zeta/\mu$','Interpreter','latex');

f = 9;style = styles{f};
load('b9.mat');
scatter(drag2fric,b,style);
legend("DP, $K_b=0.1$","DP, $K_b=100$", "SP, $K_{sp}=50$","SP, $K_{sp}=1000$","SP, $K_{sp}=20$",'Interpreter','latex');
xlim([10^-6 10^5])
ylim([0.4 1.6])
%% sp curve
clc
clear
close all
styles = {'+','o','*','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','c','m','y'};

figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


figure(4);hold on;box on;
set(gcf,'color','w');
ax = gca;
ax.XScale = "log";
ax.YScale = "log";
xlabel('$w-k$','Interpreter','latex');
ylabel('$Q/C$','Interpreter','latex');
inc = log10(1.4)-0.5*log10(2.5);
yy=10^(0.5*log10(9)+inc);
plot([2.5,9],[1.4,yy],'Color','k','LineStyle','--','HandleVisibility','off')

i = 1;
basefolder = "~/project/cells321/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(1,8,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells318/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,19,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells320/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,19,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells327/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,19,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells329/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,19,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);

set(gca,'FontSize',15);
for j = 3:4
figure(j)
legend("$k=10$","$k=100$","$k=1e3$","$k=1e4$","$k=1e5$",'Interpreter','latex');
set(gca,'FontSize',15);
end

figure(3)
xlim([1 12])
ylim([0.008 0.03])

figure(4)
xlim([1 10])
ylim([1 10])

%% dpm curve
clc
clear
close all
styles = {'+','o','*','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','c','m','y'};

figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


figure(4);hold on;box on;
set(gcf,'color','w');
ax = gca;
ax.XScale = "log";
ax.YScale = "log";
xlabel('$w-k$','Interpreter','latex');
ylabel('$Q/C$','Interpreter','latex');
inc = log10(1.4)-0.5*log10(2.5);
yy=10^(0.5*log10(9)+inc);
plot([2.5,9],[1.4,yy],'Color','k','LineStyle','--','HandleVisibility','off')
% inc = log10(3)-1.5*log10(2.5);
% yy=10^(1.5*log10(9)+inc);
% plot([2.5,9],[3,yy],'Color','k','LineStyle','--')

i = 1;
basefolder = "~/project/cells302/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells312/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells313/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,10,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells314/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,19,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);

for j = 3:4
figure(j)
legend("$K_b=0.1$","$K_b=1$","$K_b=10$","$K_b=100$","$k=1e5$",'Interpreter','latex');
set(gca,'FontSize',15);
end

figure(4)
xlim([1 10])
ylim([1 10])

figure(3)
xlim([1 12])
ylim([0.008 0.03])

%% 1/2 curve

clc
clear
close all
styles = {'+','o','*','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','c','m','y'};

figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w/\sigma_{avg}$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


figure(4);hold on;box on;
set(gcf,'color','w');
ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
xlabel('$\log_{10}(w/\sigma_{avg}-k)$','Interpreter','latex');
ylabel('$\log_{10}(Q/C)$','Interpreter','latex');
inc = log10(1.1)-0.5*log10(1.5);
yy=10^(0.5*log10(9)+inc);
plot(log10([1.5,9]),log10([1.1,yy]),'Color','k','LineStyle','--','HandleVisibility','off')
% inc = log10(3)-1.5*log10(2.5);
% yy=10^(1.5*log10(9)+inc);
% plot([2.5,9],[3,yy],'Color','k','LineStyle','--')

i = 1;
basefolder = "~/project/cells302/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(1,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells313/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,10,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells318/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(3,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells327/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,19,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells317/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(3.5,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

for j = 3:4
figure(j)
% legend("frictionless DP, $K_b=10^{-1}$","frictionless DP, $K_b=10$", "SP, $K_{sp}=10^2$","SP, $K_{sp}=10^4$", "bumpy DP",'Interpreter','latex');
set(gca,'FontSize',15);
end

% figure(4)
% xlim([1 10])
% ylim([1 10])

figure(3)
xlim([2 12])
% ylim([0.008 0.03])
%% 3/2 curve

clc
clear
close all
styles = {'+','o','*','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','c','m','y'};

figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w/\sigma_{avg}$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


figure(4);hold on;box on;
set(gcf,'color','w');
ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
xlabel('$\log_{10}(w/\sigma_{avg}-k)$','Interpreter','latex');
ylabel('$\log_{10}(Q/C)$','Interpreter','latex');
inc = log10(3)-1.5*log10(2.5);
yy=10^(1.5*log10(9)+inc);
plot(log10([2.5,9]),log10([3,yy]),'Color','k','LineStyle','--','HandleVisibility','off')

i = 1;
basefolder = "~/project/cells294/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(1,18,2,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells361/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,25,15);
cli.HopperProperty(18, 99, 18, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(1,18,2,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells311/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(3,18,2,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells362/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(18, 99, 18, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(3,18,2,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells316/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(3,18,2,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

for j = 3:4
figure(j)
% legend("frictionless DP, $K_b=10^{-1}$","frictionless DP, $K_b=10$", "SP, $K_{sp}=10^2$","SP, $K_{sp}=10^4$", "bumpy DP",'Interpreter','latex');
set(gca,'FontSize',15);
end
figure(4)
xlim([0.2 1.1])
% ylim([1 40])
%% continious flow
clc; close all; clear;
basefolder = "~/project/cells330/";
cli = CLI_hopper(basefolder);
cli.HopperProperty(0, 9, 0, 9, 3, 50);
flow = cli.hopperProperty(1).flow * 10 * 0.1 * (pi / 4) * sqrt(1/0.1)/0.005;
meanQ = mean(flow(floor(length(flow)/10):end),'all');
t = (1:length(flow))*1e4*0.005/sqrt(10);
figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$t$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');
plot(t,flow);
yline(meanQ,'LineStyle','--')
set(gca,'FontSize',15);
%% rigid limit

clc
clear
close all
styles = {'+','*','o','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','m','c','y'};

figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w/\sigma_{avg}$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


figure(4);hold on;box on;
set(gcf,'color','w');
ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
xlabel('$log(w-k)$','Interpreter','latex');
ylabel('$log(Q/C)$','Interpreter','latex');
inc = log10(1.4)-0.5*log10(2.5);
yy=10^(0.5*log10(9)+inc);
plot(log10([2.5,9]),log10([1.4,yy]),'Color','k','LineStyle','--','HandleVisibility','off')
% inc = log10(3)-1.5*log10(2.5);
% yy=10^(1.5*log10(9)+inc);
% plot([2.5,9],[3,yy],'Color','k','LineStyle','--')

i = 1;
basefolder = "~/project/cells302/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(1,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells314/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(5,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;


basefolder = "~/project/cells318/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(1,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

basefolder = "~/project/cells327/";
cli = CLI_hopper(basefolder);cli.setInterval(0.5,0.1,35,15);
cli.HopperProperty(0, 99, 0, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(5,18,1,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

for j = 3:4
figure(j)
% legend("frictionless DP, $K_b=0.1$","frictionless DP, $K_b=100$","SP, $K_{sp}=1e2$", "SP, $K_{sp}=1e5$",'Interpreter','latex');
set(gca,'FontSize',15);
end
figure(3)
xlim([2 11])

%% 3D 
clc
clear
close all
styles = {'+','o','*','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','c','m','y'};

figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w/\sigma_{avg}$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


figure(4);hold on;box on;
set(gcf,'color','w');
ax = gca;
% ax.XScale = "log";
% ax.YScale = "log";
xlabel('$\log_{10}(w/\sigma_{avg}-k)$','Interpreter','latex');
ylabel('$\log_{10}(Q/C)$','Interpreter','latex');
inc = log10(6)-2.5*log10(2.5);
yy=10^(2.5*log10(9)+inc);
plot(log10([2.5,9]),log10([6,yy]),'Color','k','LineStyle','--','HandleVisibility','off')
inc = log10(3)-1.5*log10(2.5);
yy=10^(1.5*log10(9)+inc);
plot(log10([2.5,9]),log10([3,yy]),'Color','k','LineStyle','--','HandleVisibility','off')


i = 1;
basefolder = "~/project/cells364/";
cli = CLI_hopper(basefolder);cli.setInterval(3,0.4,10,5);
cli.HopperProperty(17, 99, 17, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,18,2,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;  

basefolder = "~/project/cells366/";
cli = CLI_hopper(basefolder);cli.setInterval(3,0.4,15,5);
cli.HopperProperty(18, 99, 18, 0, 3, 50);
[fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,10.5,4,0);
plotFitted(fitted, width_b, flowRate_b, color, styles, i);
i = i+1;

% basefolder = "~/project/cells366/";
% cli = CLI_hopper(basefolder);cli.setInterval(3,0.4,10,5);
% cli.HopperProperty(17, 99, 17, 0, 3, 50);
% [fitted, width_b, flowRate_b] = cli.helper_plotFlowRate1(2,10,2,0);
% plotFitted(fitted, width_b, flowRate_b, color, styles, i);
% i = i+1;  
for j = 3:4
figure(j)
% legend("3D SP, $\mu = 0,\zeta=0.1$","3D SP, $\zeta = 0,\mu=1$",'Interpreter','latex');
set(gca,'FontSize',15);
end

figure(4)
xlim([0.3 1.1])
ylim([0.2 2.7])

%% exp data
clc
clear
close all
styles = {'+','o','*','s','^','p','|','x','>','<','v','_','d'};
color = {'r','g','b','c','m','y'};

expFlow = readtable("fig5.xlsx");
dataX = table2array(expFlow(:,1));
dataX = dataX(dataX>0);
dataY = table2array(expFlow(:,2));
dataY = dataY(dataY>0);
figure(4);hold on;box on;
xlim([-0.5 1.2])
set(gcf,'color','w');
ax = gca;
inc = log10(0.7)-0.5*log10(0.5);
yy=10^(0.5*log10(9)+inc);
plot(log10([0.5,9]),log10([0.7,yy]),'Color','k','LineStyle','--','HandleVisibility','off')
xlabel('$\log_{10}(w/\sigma_{avg}-k)$','Interpreter','latex');
ylabel('$\log_{10}(Q/C)$','Interpreter','latex');
figure(3);hold on;box on;
set(gcf,'color','w');
xlabel('$w/\sigma_{avg}$','Interpreter','latex');
ylabel('$Q$','Interpreter','latex');


% factor = 1;

b_fit = 0.5;
i=1;
fitfun = fittype( @(a, c, x) abs(c*(x - a).^b_fit));
fitted = fit( dataX, dataY, fitfun, 'StartPoint', [0.01,1.1]);

g_eff = 9.8 * (997-963)/963;
factor = (pi / 4) * sqrt(3e-4/g_eff);
fitted.c = fitted.c * factor;
plotFitted(fitted, dataX, factor*dataY, color, styles, i);
figure(3)
xlim([1 14])
% ylim([0.2 2.7])

for j = 3:4
figure(j)
% legend("3D SP, $\mu = 0,\zeta=0.1$","3D SP, $\zeta = 0,\mu=1$",'Interpreter','latex');
set(gca,'FontSize',15);
end

%%
function plotFitted(fitted, width_b, flowRate_b, color, styles, i)

factor = 0.1 * (pi / 4) * sqrt(1/0.1)/0.005;
factor = 1;
fitted.c = fitted.c * factor;
flowRate_b = flowRate_b * factor;
fittedX = linspace(min(width_b,[],'all'),max(width_b,[],'all'),100);
fittedY = feval(fitted, fittedX);
figure(3);plot(fittedX, fittedY,'k','LineWidth',1,'HandleVisibility','off');
scatter(width_b, flowRate_b,append(color{i},styles{i}));
figure(4);scatter(log10(width_b-fitted.a),log10(abs(flowRate_b/fitted.c)),append(color{i},styles{i}));
end