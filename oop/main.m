clc
clear
close all
addpath FlowRate
addpath Trial
addpath CLI

basefolder = "~/project/cells172/";
%basefolder = "~/project/cells142/";
% 
% %folderList = ["~/project/cells102/","~/project/cells101/","~/project/cells100/","~/project/cells94/","~/project/cells98/"];
% %folderList = ["~/project/cells105/","~/project/cells104/","~/project/cells103/","~/project/cells106/"];
% % % 
cli = CLI_DPM(basefolder);
%cli.readSysProperty(0, 5, 0, 5);
cli.readTao();

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

cli = CLI_hopper(basefolder);
%cli.calHopperProperty(99, 39, 4, 50);
% cli.plotClogP();
% folderList = ["~/project/cells145/","~/project/cells146/","~/project/cells143/","~/project/cells144/","~/project/cells133/","~/project/cells123/"];
% cli.compare(folderList);
% cli.calHopperProperty(39, 19, 1);
% cli.plotFlowRate();
% cli.calHopperProperty(39, 19, 2, 10);
% cli.plotFlowRate();
% cli.calHopperProperty(39, 19, 2, 50);
% cli.plotFlowRate();

cli.HopperProperty(99, 39, 0, 0, 3, 50);
% cli.plotFlowRateWithN()
% cli.plotFlowRateVSN_W()