function [result] = TemperatureCalculation(datafilepath, calfilepath, measfilepath, measfilename, ...
    pathlength, pressure, yh2o, i, leftStart, rightEnd, broadbandType)

% load conversion coefficients (variables: 'mus', 'coeffs')
load(strcat(calfilepath,'mvToCvCoeffs.mat'));
%load('F:\John Tobiasson\MATLAB Files\johnsAttempt\currentCode\mvToCvCoeffs-02-17-2016-1.mat');

%% Open file, identify wavenumber band locations, prune measured values to size,
% convert to intensities

% Get the raw test data
% if i<10 %matches file format "...0001.csv" from OMNIC default
%     numstr = strcat('000',num2str(i));
% elseif i<100
%     numstr = strcat('00',num2str(i));
% elseif i<1000
%     numstr = strcat('0',num2str(i));
% else
numstr = num2str(i);
% end

rawData = csvread(strcat(measfilepath,measfilename,numstr,'.csv'));
rawData(:,1) = rawData(:,1)+(rawData(:,1)*3E-5 - .0143); %measured to modeld offset from peak analysis
% + .167; %JT shift. Wrong, but better than nothing

% Get bounding rows for left and right BG and A-E ratio wavelength bands
rawWavenumbers = rawData(:,1);
[~, removeLow] = min(abs(rawWavenumbers-leftStart));
[~, removeHigh] = min(abs(rawWavenumbers-rightEnd));

%Cut off unnecessary data and re-assign matrices
nuMeas = rawData(removeLow:removeHigh, 1);
voltages = rawData(removeLow:removeHigh, 2);

% Get data bounds
[~, rowLowA] = min(abs(nuMeas-5185));
[~, rowLowB] = min(abs(nuMeas-5310));
[~, rowLowC] = min(abs(nuMeas-5435));
[~, rowHighC] = min(abs(nuMeas-5560));
[~, rowLowE] = min(abs(nuMeas-5615));
[~, rowHighE] = min(abs(nuMeas-5715));

[~, rowLowLeft] = min(abs(nuMeas-(4630-25))); %4600-4750 old BG
[~, rowHighLeft] = min(abs(nuMeas-(4630+25)));

[~, rowLowRight] = min(abs(nuMeas-(6150-25))); %5800-5900 old BG
[~, rowHighRight] = min(abs(nuMeas-(6150+25)));

% Compile band wavenumber indexes into array for BG calc
bandsBroadband = [rowLowLeft, rowHighLeft, rowLowRight, rowHighRight];

% convert all Mv's directly to I's (W/m^2/sr/um)
Cvs = MvToCv(nuMeas, voltages, coeffs, mus);
intsMeas = (voltages ./ Cvs);

%% Calculate background radiation

h = 6.62606957 * (10^(-34));
k = 1.3806488 * (10^(-23));
c0 = 2.99792458 * (10^8);

[tempBroadband, emBroadband, ~, ~] ...
    = BackgroundTempCalculation(bandsBroadband, nuMeas, intsMeas);

tempGasPrevious = 1300; %initial gas temperature guess
tempGasAvg = tempGasPrevious;
intsMeasNoBG = ones(length(nuMeas), 1); %for first iteration
tempChange = 1; %change between current and previous temperature calculations
iter = 0; %iteration counter

%% Get the gas temperature
while abs(tempChange) >= .0001 && iter <= 50
    
    % calculate gas intensities without background emission. Iterates for absorption coefficients
    intsMeasBG = (2 * h * (c0^2) .* (nuMeas.^5) * (10^4)) ...
        ./ ( exp((h * c0 .* nuMeas * 100) / (k .* tempBroadband)) - 1 ) ...
        * emBroadband; %broadband emission intensity (particle or wall)
    
    gasBBIntensity = (2 * h * (c0^2) .* (nuMeas.^5) * (10^4)) ...
        ./ ( exp((h * c0 .* nuMeas * 100) / (k * tempGasAvg)) - 1 ); %BB gas emission
    
    gasEmissivity = intsMeasNoBG ./ gasBBIntensity; %spectral measured emissivity
    
    if strcmp(broadbandType, 'wall') == 1 || strcmp(broadbandType, 'removeWall') == 1
        
        absorptionCoeffs = 1 - gasEmissivity;
        
    elseif strcmp(broadbandType, 'particle') == 1
        
        absorptionCoeffs = 1; % assume no interaction of particles and gas
        
    else % unknown broadband radiation type
        
        broadbandType = throwerror; %this will throw an error if executed
        
    end
    
    intsMeasNoBG = intsMeas - (intsMeasBG .* absorptionCoeffs);
    
    % get intensities for each band
    intsMeasA = intsMeasNoBG(rowLowA:rowLowB);
    intsMeasB = intsMeasNoBG(rowLowB:rowLowC);
    intsMeasC = intsMeasNoBG(rowLowC:rowHighC);
    intsMeasE = intsMeasNoBG(rowLowE:rowHighE);
    
    % get wavelength band values
    nuMeasA = 10^4 ./ nuMeas(rowLowA:rowLowB);
    nuMeasB = 10^4 ./ nuMeas(rowLowB:rowLowC);
    nuMeasC = 10^4 ./ nuMeas(rowLowC:rowHighC);
    nuMeasE = 10^4 ./ nuMeas(rowLowE:rowHighE);
    
    % Change sums to integrals, flipping arrays so the integrals are +
    intMeasA = trapz(flipud(nuMeasA), flipud(intsMeasA));
    intMeasB = trapz(flipud(nuMeasB), flipud(intsMeasB));
    intMeasC = trapz(flipud(nuMeasC), flipud(intsMeasC));
    intMeasE = trapz(flipud(nuMeasE), flipud(intsMeasE));
    
    % calcuate ratio bands for temperature curve fit
    ratioMeasEA = intMeasE / intMeasA;
    ratioMeasEB = intMeasE / intMeasB;
    ratioMeasEC = intMeasE / intMeasC;
    
    % Get temperature from each band
    load('C:\Users\Scott\Documents\MATLAB\CabsDatabaseFiles\ratioModel') %all ratios
    % variables 'ratioModelEA','ratioModelEB','ratioModelEC'
    
    pathlengthModel = (.02:.02:1)'; %data points present in the files (same order as indexed in array)
    pressureModel = [.5, 1, 2, 4, 8, 15, 30, 40]';
    yh2oModel = [.05, .10, .20, .30, .40]';
    tempModel = (300:100:3000);
    
    for i = 1:length(tempModel) %interpolate known conditions to find ratioModelEA @ conditions = f(Temperature)
        
        ratioModelEAofT(i) = interpn(pathlengthModel, pressureModel, yh2oModel,...
            ratioModelEA(:,:,:,i,:), pathlength, pressure, yh2o, 'spline');
        ratioModelEBofT(i) = interpn(pathlengthModel, pressureModel, yh2oModel,...
            ratioModelEB(:,:,:,i,:), pathlength, pressure, yh2o, 'spline');
        ratioModelECofT(i) = interpn(pathlengthModel, pressureModel, yh2oModel,...
            ratioModelEC(:,:,:,i,:), pathlength, pressure, yh2o, 'spline');
        
    end
    
    %interpolate temperature = f(ratioModelEA @ conditions) to find temperature @ EA measured
    tempEA = interp1(ratioModelEAofT, tempModel, ratioMeasEA, 'spline');
    tempEB = interp1(ratioModelEBofT, tempModel, ratioMeasEB, 'spline');
    tempEC = interp1(ratioModelECofT, tempModel, ratioMeasEC, 'spline');
    
    % Average temperatures
    tempGasAvg = mean([tempEC, tempEB, tempEA]);
    
    if abs(1400 - tempGasAvg) > 500 %verify temperature is reasonable (define reasonable...)
        
        tempGasAvg = throwerror;
    end
    
    tempChange = (tempGasAvg - tempGasPrevious) / tempGasPrevious;
    tempGasPrevious = tempGasAvg; % update 'previous' value for next iteration
    iter = iter + 1;
end

% Get the signal to noise ratio
%[~, waterToBB] = NoiseAnalysis(rawData, intsMeas, backgroundEmissivity, backgroundTemp, nuMeas);

%Calculate average gas emissivity in H2O region
eEmis = mean(gasEmissivity(rowLowE:rowHighE));

%% Compile and save results
result(1) = tempEA;
result(2) = tempEB;
result(3) = tempEC;
result(4) = tempGasAvg;
result(5) = tempBroadband;
result(6) = emBroadband;
result(7) = eEmis;

% close all
% plot(nuMeas,[intsMeas, intsMeasNoBG, intsMeasBG])
% ylabel('emission intensity [W/m^2/sr/um]')
% xlabel('wavenumber [cm-1]')
% legend('intsMeas', 'intsMeasNoBG', 'intsMeasBG')

% save the measured data
save(strcat(measfilepath,measfilename,numstr,'meas.mat'), 'nuMeas', 'intsMeas', 'intsMeasNoBG', 'intsMeasBG')

pause

end
