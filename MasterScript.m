%% Setup
close all
clear variables
clc
pause off


%% System Calibration
calfilepath = 'C:\Users\Scott\Documents\Working Files\All BFR Combustion\extended spectrum\April Calibration\';
calfilename = '04282016Callibration'; %code assumes "...0001.csv" tag on end of file
filterfilename = '04282016CalibrationFiltered';
datafilepath = 'F:\';

starttemp = 600; %max BB calibration temp
endtemp = 1100; %min BB calibration temp
tempstep = 50; %delta temperature 

leftStart = 3750; %4600-5900 range for T with BG
rightEnd = 7800;  %4500-7500 good calibration range
                  %3750-7800 good filter range

%FilterData(calfilepath, calfilename, filterfilename, datafilepath,...
%    starttemp, endtemp, tempstep, leftStart, rightEnd);

leftStart = 4500; %4600-5900 range for T with BG
rightEnd = 7500;  %4500-7500 good calibration range
                  %3750-7800 good filter range

%CalibrationCurves(calfilepath, filterfilename, tempStart, tempEnd, tempStep, leftStart, rightEnd)
%PlotCalibratedData(calfilepath, filterfilename, starttemp, endtemp, tempstep, leftStart, rightEnd)

%% Data Processing

measfilepath = 'C:\Users\Scott\Documents\Working Files\All BFR Combustion\extended spectrum\IFFW\';
date = '04-29-2016';
measfilename = strcat(date,'-Test');
opticalpathfilename = strcat(date,'-Path');
suctionpyrofilename = strcat(date,'-SuctionPyro');

broadbandType = 'particle'; % 'wall' 'removeWall' 'particle'
fuelType = 'fine wood'; %'nat gas' 'fine wood' 'med wood' 

startdatafile = 2;
stopdatafile = 50;
stepdatafile = 1;

pressure = .844;

MasterProcessData(datafilepath, calfilepath, measfilepath, measfilename, opticalpathfilename,...
   startdatafile, stopdatafile, stepdatafile, pressure, leftStart, rightEnd, broadbandType, fuelType)

%MasterSpectralData(calfilepath, measfilepath, measfilename, suctionpyrofilename,...
%    startdatafile, stopdatafile, stepdatafile, pressure, leftStart, rightEnd)
%CompareModels(pressure, leftStart, rightEnd) 



%% ----------------Documentation-------------------
%{

% For the 04282016 calibration (PFFW, IFFW)

calfilepath = 'C:\Users\Scott\Documents\Working Files\All BFR Combustion\extended spectrum\April Calibration\';
calfilename = '04282016Callibration'; %code assumes "...0001.csv" tag on end of file
filterfilename = '04282016CalibrationFiltered';
datafilepath = 'F:\';

starttemp = 600; %max BB calibration temp
endtemp = 1100; %min BB calibration temp
tempstep = 50; %delta temperature 

% For the 02172016 calibration (PFNG, IFNG, PFMW)

calfilepath = 'C:\Users\Scott\Documents\Working Files\All BFR Combustion\extended spectrum\February Calibration\';
calfilename = '02172016Callibration'; %code assumes "...0001.csv" tag on end of file
filterfilename = '02172016CalibrationFiltered';
datafilepath = 'F:\';

starttemp = 1100; %max BB calibration temp
endtemp = 600; %min BB calibration temp
tempstep = -50; %delta temperature 

%}

