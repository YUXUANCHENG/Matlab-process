clc
clear
close all
addpath FlowRate

basefolder = "~/project/cells97/";
%basefolder = "~/project/cells102/";

folderList = ["~/project/cells102/","~/project/cells94/","~/project/cells98/"];
%folderList = ["~/project/cells105/","~/project/cells104/","~/project/cells103/"];
CLI.compare_DPM_SP(folderList);

% cli = CLI_hopper(basefolder);
% %cli.calHopperProperty(1, 6, 3, 10);
% % cli.calHopperProperty(39, 19, 1);
% % cli.plotFlowRate();
% 
% cli.HopperProperty(39, 19, 0, 19, 3, 10);
% cli.plotFlowRateWithN()