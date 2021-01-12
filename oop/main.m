clc
clear
close all
addpath FlowRate
basefolder = "~/project/cells97/";
% fileReader = FileReader();
% trial = Trial_hopper(1,3,basefolder,fileReader);
% %trial.plotInitial();
% trial.readMDdata();
% trial.plotLastFrame(2);
% trial.createCalculator();
% trial.printCellCount();
% trial.flowRate(2,10);

cli = CLI_hopper(basefolder);
%cli.hopperProperty(39, 19);
%cli.calHopperProperty(2, 5, 2, 10);
cli.calHopperProperty(39, 19, 1);
cli.plotFlowRate();