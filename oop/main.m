clc
clear
close all
addpath FlowRate
%basefolder = "~/project/cells97/";
basefolder = "~/project/cells102/";

fileReader = FileReader_back();
trial = Trial(1,3,basefolder,fileReader);
%trial.plotInitial();
trial.readMDdata();
trial.plotLastFrame(2);
trial.createCalculator();
trial.plotVelDistribution();
trial.cal_msd();
trial.plotMSD();

% cli = CLI_hopper(basefolder);
% %cli.hopperProperty(39, 19);
% %cli.calHopperProperty(2, 5, 2, 10);
% cli.calHopperProperty(39, 19, 1);
% cli.plotFlowRate();