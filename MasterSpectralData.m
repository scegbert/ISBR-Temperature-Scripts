function [] = MasterSpectralData(calfilepath, measfilepath, measfilename, suctionpyrofilename,...
    startdatafile, stopdatafile, stepdatafile, pressure, leftStart, rightEnd)

%compares measured intentities, models created from optically measured intensities,
%and models created from suction pyrometer values.

leftStartModel = 200;     %expanded range to obtain the full model
rightEndModel = 9600;  %max range for Pearson data set is 0-10,000 cm-1

%read in optical data results
results = xlsread(strcat(measfilepath,measfilename,'-Results'));

TGasOptical = results(:,7);
backgroundTemp = results(:,8); %ones(length(TGasOptical),1); <----- change for BG conditions
backgroundEmissivity = results(:,9); %zeros(length(TGasOptical),1); <----- change for BG conditions
pathLengths = results(:,2);
YH2Oopt = results(:,11);
YCO2opt = results(:,12);

%read in suction pyro and exhaust gas measurements
suctionpyro = xlsread(strcat(measfilepath,suctionpyrofilename,'.xlsx'));

TSuctionPyro = suctionpyro(:,2);
YH2Oexh = suctionpyro(:,3);

HtoC = 4; % natural gas
H2OtoC = 0;
YCO2exh = YH2Oexh / (0.5 * HtoC + H2OtoC);

% load conversion coefficients (measured values to intensities)
load(strcat(calfilepath,'mvToCvCoeffs.mat'));

sensT = 300;
sensEm = 1;
partT = 1;
partKappa = 0;

bands = [200, 1000;
    1000, 2100;
    2100, 2600;
    2600, 4600;
    4600, 6000;
    6000, 9600; %major bands
    4700, 4850;
    4850, 5025;
    5025, 5185; 
    5185, 5310;
    5310, 5435;
    5435, 5560; %Broad - A
    5615, 5715;
    5715, 5850; %E - Broad
    6800, 7000]; %Fringe

% [4750, 5020;
%     5020, 5185;
%     5185, 5310; %A
%     5310, 5435; %B
%     5435, 5560; %C
%     5615, 5715; %E
%     6600, 6700;
%     6700, 6800;
%     6800, 6900;
%     6900, 7000;
%     7000, 7100;
%     7100, 7200;
%     7200, 7300;
%     4600, 5900; %temperature band
%     leftStart, rightEnd; %all measured
%     leftStartModel, rightEndModel]; %all modeled

for i = startdatafile:stepdatafile:stopdatafile
    
    i
    tic
    typestr = 'opt';
    
    intsAllOpt = GenerateSpectra(measfilepath, measfilename, i, typestr, TGasOptical(i),...
        YH2Oopt(i), YCO2opt(i), partT, partKappa, backgroundTemp(i),...
        backgroundEmissivity(i), sensT, sensEm, pressure, pathLengths(i),...
        leftStartModel, rightEndModel, bands);
    
    typestr = 'spe';
    
    intsAllSPE = GenerateSpectra(measfilepath, measfilename, i, typestr, TSuctionPyro(i),...
        YH2Oexh(i), YCO2exh(i), partT, partKappa, backgroundTemp(i),...
        backgroundEmissivity(i), sensT, sensEm, pressure, pathLengths(i),...
        leftStartModel, rightEndModel, bands);
    
    intsAllSPE(:,1) = []; % measured values are already in the other matrix
    
    intsAll = horzcat(bands, intsAllOpt, intsAllSPE);
    
    headers = {'start cm-1', 'stop cm-1',...
        'Measured', 'Modeled Opt', 'Convolved Opt',...
        'Modeled H2O only Opt', 'Convolved H2O only Opt',...
        'Modeled CO2 only Opt', 'Convolved CO2 only Opt',...
        'Modeled SP&E', 'Convolved SP&E',...
        'Modeled H2O only SP&E', 'Convolved H2O only SP&E',...
        'Modeled CO2 only SP&E', 'Convolved CO2 only SP&E'};
    
    numstr = num2str(i);
    
    xlswrite(strcat(measfilepath,measfilename,numstr,'-Integrated'), headers, 1, 'A1');
    xlswrite(strcat(measfilepath,measfilename,numstr,'-Integrated'), intsAll, 1, 'A2');
    
end

end

