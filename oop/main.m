clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI

%basefolder = "~/project/cells95/";
basefolder = "~/project/cells133/";

%folderList = ["~/project/cells102/","~/project/cells101/","~/project/cells100/","~/project/cells94/","~/project/cells98/"];
%folderList = ["~/project/cells105/","~/project/cells104/","~/project/cells103/","~/project/cells106/"];
% 
cli = CLI_DPM([]);
%folderList = ["~/project/cells134/","~/project/cells135/","~/project/cells113/"];
folderList = ["~/project/cells107/"];
%folderList = ["~/project/cells128/","~/project/cells127/", "~/project/cells130/", "~/project/cells129/"];
%folderList = ["~/project/cells113/","~/project/cells107/", "~/project/cells110/", "~/project/cells111/"];
cli.compare(folderList);
% % cli = CLI_Disk([]);
% % folderList = ["~/project/cells116/", "~/project/cells110/"];
% % %folderList = ["~/project/cells117/"];
% % cli.compare(folderList);

% cli = CLI_hopper(basefolder);
% cli.calHopperProperty(99, 39, 4, 50);
% cli.plotClogP();
% cli.calHopperProperty(39, 19, 1);
% cli.plotFlowRate();
% cli.calHopperProperty(39, 19, 2, 10);
% cli.plotFlowRate();
% cli.calHopperProperty(39, 19, 2, 50);
% cli.plotFlowRate();

%cli.HopperProperty(99, 39, 0, 38, 4, 50);
% cli.plotFlowRateWithN()
% cli.plotFlowRateVSN_W()