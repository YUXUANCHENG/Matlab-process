clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI

% % basefolder = "~/project/cells245/";
% basefolder = "~/scratch60/cells245/";
% % 
% % %folderList = ["~/project/cells102/","~/project/cells101/","~/project/cells100/","~/project/cells94/","~/project/cells98/"];
% % %folderList = ["~/project/cells105/","~/project/cells104/","~/project/cells103/","~/project/cells106/"];
% % % % 
% cli = CLI_DPM(basefolder);
% cli.readSysProperty(10, 0, 9, 0);
% %cli.Arrhenius(0)
% %cli.readTao(2, 0.787, 1.4, 1.8, 0.1)
% cli.plotPhaseDiagram();
% cli.Angell(1);

% folderList = ["~/project/cells179/","~/project/cells180/","~/project/cells184/","~/project/cells182/","~/project/cells183/"];
% for basefolder = folderList
%     cli = CLI_DPM(basefolder);
%     cli.Angell(1);
% end

% cli.plotScalling();
%folderList = ["~/project/cells134/","~/project/cells135/","~/project/cells113/"];
%folderList = ["~/project/cells107/"];
%folderList = ["~/project/cells139/","~/project/cells137/","~/project/cells129/"];
%folderList = ["~/project/cells76/","~/project/cells98/", "~/project/cells75/"];
%folderList = ["~/project/cells113/","~/project/cells107/", "~/project/cells110/", "~/project/cells111/"];
%cli.compare(folderList);
% % cli = CLI_Disk([]);
% % folderList = ["~/project/cells116/", "~/project/cells110/"];
% % %folderList = ["~/project/cells117/"];
% % cli.compare(folderList);
basefolder = "~/project/cells342/";
%basefolder = "~/project/test3/";
cli = CLI_hopper(basefolder);
% cli.HopperProperty(0, 99, 0, 69, 3, 50);
% cli.HopperProperty(0, 14, 0, 3, 3, 50);
% % % % cli.calHopperProperty(29, 39, 3, 50);
% % % %cli.calHopperProperty(0, 10, 3, 50);
% cli.plotClogP();

%folderList = ["~/project/cells145/","~/project/cells146/","~/project/cells143/","~/project/cells144/","~/project/cells133/","~/project/cells123/"];
%folderList = ["~/project/cells154/","~/project/cells156/","~/project/cells155/","~/project/cells157/"];
% folderList = ["~/project/cells95/","~/project/cells97/"];
% cli.compare(folderList);
% cli.calHopperProperty(39, 19, 1);
% cli.plotFlowRate();
% cli.calHopperProperty(39, 19, 2, 10);
% cli.plotFlowRate();
% cli.calHopperProperty(39, 19, 2, 50);
% cli.plotFlowRate();

cli.HopperProperty(16, 99, 16, 10, 0, 50);
% cli.calHopperProperty(0, 99, 1);
cli.helper_plotFlowRate1(3,20,2,1);
% cli.plotFlowRateWithN()
% cli.plotFlowRateVSN_W()