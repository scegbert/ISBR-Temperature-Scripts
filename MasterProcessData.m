function [results] = MasterProcessData(datafilepath, calfilepath, measfilepath, measfilename, opticalpathfilename,...
    startdatafile, stopdatafile, stepdatafile, pressure, leftStart, rightEnd, broadbandType, fuelType)

%read in optical pathlength data
opticalpathfile = xlsread(strcat(measfilepath,opticalpathfilename));

opticalPos = .01 * opticalpathfile(:,2); %position of optical probe
BBPos = .01 * opticalpathfile(:,3); %position of BB target in BFR

pathLengths = 0.75 - (opticalPos + BBPos); % This is for BFR data at BYU
pathCenters = (pathLengths / 2) + opticalPos;

H2OWavenumbers = [5615, 5715]; %this is band E, compareable to bands A-C [5185 5560];
yh2o = .16; %initial guess at yh20

for i = startdatafile:stepdatafile:stopdatafile
    
    i
    tic
        
    resultsI(1:3) = [i, pathLengths(i), pathCenters(i)];
    
    yh2oChange = 1; %change of Yh2o between iterations of temperature and yh2o calculations
    
    while abs(yh2oChange) > .01 %iterate until changes in gas temp don't change Yh2o
        
        %% calculate temperature
        resultsT(1:7) = TemperatureCalculation(datafilepath, calfilepath, measfilepath, measfilename, ...
            pathLengths(i), pressure, yh2o, i, leftStart, rightEnd, broadbandType);
        
        %% organize 1D RTE input conditions dependant on combustion and fuel type
        
        tempGasAvg = resultsT(4);
        
        numstr = num2str(i);
        load(strcat(measfilepath,measfilename,numstr,'meas.mat'));
        % variables: 'nuMeas', 'intsMeas', 'intsMeasNoBG', 'intsMeasBG'
        
        %determine wall and particle conditions for 1D RTE
        if strcmp(broadbandType, 'wall') == 1
            
            wallConditions = [resultsT(5),resultsT(6)]; % use calculated Wall T, Wall em
            particleConditions = [300, 0]; % ignore Particle T, Particle Kappa
            
            intMeasH2O = IntegrateWavelength(intsMeas, nuMeas, H2OWavenumbers(1), H2OWavenumbers(2));
            
        elseif strcmp(broadbandType, 'removeWall') == 1
            
            wallConditions = [300, 0]; % ignore Wall T, Wall em
            particleConditions = [300, 0]; % ignore Particle T, Particle Kappa
            
            intMeasH2O = IntegrateWavelength(intsMeasNoBG, nuMeas, H2OWavenumbers(1), H2OWavenumbers(2));
            
        elseif strcmp(broadbandType, 'particle') == 1
            
            wallConditions = [300, 0]; % ignore Wall T, Wall em
            
            partKappa = -log(1-resultsT(6)) / pathLengths(i); % results(i,10) = partEmissivity
            particleConditions = [resultsT(5), partKappa]; %use calculated Particle T, Particle Kappa
            
            intMeasH2O = IntegrateWavelength(intsMeas, nuMeas, H2OWavenumbers(1), H2OWavenumbers(2));
            
        else % unknown broadband radiation type
            broadbandType = throwerror; %unknown broadband type
            
        end
        %% calculate Y H2O
        
        [resultsY(1), resultsY(2)] = H2OCalculation(intMeasH2O, tempGasAvg, wallConditions,...
            particleConditions, pressure, yh2o, pathLengths(i), H2OWavenumbers, fuelType);
        
        yh2oChange = (resultsY(1) - yh2o) / yh2o;
        yh2o = resultsY(1); %update yh2o value used in temperature calculations
            %note that 'resultsY(1)' NOT 'yh2o' is recorded as the final concentration value
            %these two values should be nearly identical before exiting the while loop
    end
    
    ProcessDataTime = toc
    
    results(i,:) = horzcat(resultsI, resultsT, resultsY, ProcessDataTime); %combine results
    clear resultsI resultsT resultsY ProcessDataTime
    
    %% save the data
    % saves every iteration to ensure data isn't lost if an error is found
    headers = {'Test', 'Path Length [m]', 'Path Center [m]', 'Gas T E/A [K]',...
        'Gas T E/B [K]', 'Gas T E/C [K]', 'Gas T Final [K]',...
        'Broadband T [K]', 'Broadband Epsilon', 'Band E Gas Epsilon',...
        'Y H2O', 'Y CO2', 'Time to Process [sec]'};
    
    xlswrite(strcat(measfilepath,measfilename,'-Results'), headers, 1, 'A1');
    xlswrite(strcat(measfilepath,measfilename,'-Results'), results, 1, 'A2');
    
end

end
