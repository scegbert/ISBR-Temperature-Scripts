function [] = CompareModels(pressure, leftStart, rightEnd) 

%calls oneD_RTE_GAP function
%similar to SpectralData function, with a focus on models and not comparing
%to measured as well. 

gasT = [300;1400;300;1400;1400;1400;1400;300;1400;300];


gasConcH2O = .178;
gasConcCO2 = gasConcH2O/2;

wallT = 1170;
wallEm = .04;

partT = 1; %no particles
partKappa = 0;

sensT = 1; %sensor not an emitter
sensEm = 1;

pathLength = .5; %m

startWavenumber = leftStart; %4600 - 5900 cm-1 includes background
stopWavenumber = rightEnd;  %5185 - 5715 cm-1 includes A-E
wavenumbers = transpose(startWavenumber:0.005:stopWavenumber);

num = 1; %used for evaluating multiple conditions
checkModel = 1; %1 runs the first model back through the equations to get a second model

for i = 1:1:num
    i
    
    modInts(:,i) = OneD_RTE_GAP(gasT, gasConcH2O, gasConcCO2, partT, partKappa,...
        wallT, wallEm, sensT, sensEm, pressure, pathLength, ...
        startWavenumber, stopWavenumber);

    if checkModel == 1
       
        % Get data bounds
        [~, rowLowA] = min(abs(wavenumbers-5185));
        [~, rowLowB] = min(abs(wavenumbers-5310));
        [~, rowLowC] = min(abs(wavenumbers-5435));
        [~, rowHighC] = min(abs(wavenumbers-5560));
        [~, rowLowE] = min(abs(wavenumbers-5615));
        [~, rowHighE] = min(abs(wavenumbers-5715));

        [~, rowLowLeft] = min(abs(wavenumbers-4600)); %4600-4750 BG
        [~, rowHighLeft] = min(abs(wavenumbers-4750)); 

        [~, rowLowRight] = min(abs(wavenumbers-5800)); %5800-5900 BG
        [~, rowHighRight] = min(abs(wavenumbers-5900));

        % Compile band wavenumber indexes into array for Temp calc
        TBands = [rowLowA, rowLowB, rowLowC, rowHighC, rowLowE, rowHighE];   
        BgBands = [rowLowLeft, rowHighLeft, rowLowRight, rowHighRight];

        % Compile A-C band wavenumbers into array for YH2O
        H2OBands = [rowLowA, rowHighC];
        H2OWavenumbers = [5185, 5560];
    
        pseudoRawData = [wavenumbers, modInts(:,i)];
        
        NewResults(i,1:12) = TemperatureCalculation(wavenumbers, modInts(:,i), TBands, BgBands,...
            pathLength, i, pseudoRawData);
    
        NewGasTemp = NewResults(i,5);
        NewWallConditions = [NewResults(i,9),NewResults(i,10)]; %Wall T, Wall em
    
        [NewYH2O, NewYCO2] = H2OCalculation(modInts(:,i), wavenumbers, NewGasTemp, H2OBands,...
            H2OWavenumbers, NewWallConditions, pressure, pathLength);
        
        NewModInts(:,i) = OneD_RTE_GAP(NewGasTemp, NewYH2O, NewYCO2, partT, partKappa,...
            NewWallConditions(1), NewWallConditions(2), sensT, sensEm, pressure, pathLength, ...
            startWavenumber, stopWavenumber);
        
        hold on
        figure (1)
        plot(wavenumbers,modInts(:,i),'LineWidth',2)
        plot(wavenumbers,NewModInts(:,i),'LineWidth',2)
        xlabel('Wavenumber (cm^{-1})')
        ylabel('Intensity (W/m^2/sr)')
        
        errorTgas = 100*(NewGasTemp - mean(gasT)) / mean(gasT);
        errorTwall = 100*(NewWallConditions(1) - wallT) / wallT;
        errorEmwall = 100*(NewWallConditions(2) - wallEm) / wallEm;
        errorYH2O = 100*(NewYH2O - gasConcH2O) / gasConcH2O;
        errorYCO2 = 100*(NewYCO2 - gasConcCO2) / gasConcCO2;
        
        errorTable = [mean(gasT), NewGasTemp, errorTgas;...
                      wallT, NewWallConditions(1), errorTwall;...
                      wallEm, NewWallConditions(2), errorEmwall;...
                      gasConcH2O, NewYH2O, errorYH2O]
    end
end
%title('1.88 to 1.885 um, Opt 26% larger than SP, .1% larger than Matching')
%legend('Optical (1251, .153)','Suction Pyrometer (1279, .096)','Area Matching (1354, .096)')


