function [intsAll] = GenerateSpectra(measfilepath, measfilename, i, typestr, TGas, YH2O,...
    YCO2, partT, partKappa, backgroundTemp, backgroundEmissivity,...
    sensT, sensEm, pressure, pathLength, leftStartModel, rightEndModel, bands)

%% Model with BOTH H2O and CO2
[intsModel, intsConv, nuModel] = GaussianConvolution(...
    TGas, YH2O,YCO2, partT, partKappa, backgroundTemp,...
    backgroundEmissivity, sensT, sensEm, pressure,...
    pathLength, leftStartModel, rightEndModel);

% Model with ONLY H2O, NO CO2
[intsModelH, intsConvH, ~] = GaussianConvolution(...
    TGas, YH2O, 0, partT, partKappa, backgroundTemp,...
    backgroundEmissivity, sensT, sensEm, pressure,...
    pathLength, leftStartModel, rightEndModel);

% Model with ONLY CO2, NO H2O (YH2O -> 0, got buggy with = 0)
[intsModelC, intsConvC, ~] = GaussianConvolution(...
    TGas, 1E-9, YCO2, partT, partKappa, backgroundTemp,...
    backgroundEmissivity, sensT, sensEm, pressure,...
    pathLength, leftStartModel, rightEndModel);

numstr = num2str(i);

save(strcat(measfilepath,measfilename,numstr,typestr,'.mat'),...
    'nuModel', 'intsModel', 'intsConv', 'intsModelH',...
    'intsConvH', 'intsModelC', 'intsConvC')

load(strcat(measfilepath,measfilename,numstr,'meas.mat'))
% variables: 'nuMeas', 'intsMeas', 'intsMeasNoBG', 'intsMeasBG'

% include or exclude background from Measured, as applicable
if backgroundEmissivity == 0 || backgroundTemp <=300
    intsMeas = intsMeasNoBG;
end

lineWidth = 2;

for j = 1:1:length(bands)
    
    %determine location of bands
    [~, startConv] = min(abs(nuModel-bands(j,1)));
    [~, stopConv] = min(abs(nuModel-bands(j,2)));
    
    [~, startMeas] = min(abs(nuMeas-bands(j,1)));
    [~, stopMeas] = min(abs(nuMeas-bands(j,2)));
    
    %find max value on y-axis for sizing the plot
    if max(intsMeas(startMeas:stopMeas)) > max(intsConv(startConv:stopConv))
        ymax = max(intsMeas(startMeas:stopMeas)) * 1.1;
    else
        ymax = max(intsConv(startConv:stopConv)) * 1.1;
    end
            
    % PLOT measurements, all 3 convolved models
    close all
    hold on
    plot(nuMeas, intsMeas,'LineWidth',lineWidth)
    plot(nuModel,[intsConv, intsConvH, intsConvC],'LineWidth',lineWidth)
    ylabel('Emission Intensity [W/m^2/sr/um]')
    xlabel('Wavenumber [cm-1]')
    legend('Measured Spectra', 'Model H2O and CO2',...
        'Model H2O Only', 'Model CO2 Only')
    % legend('Location','northwest')
    % ymax = ymax * 1.2;
    axis([bands(j,1) bands(j,2) 0 ymax])
    
    %saveas(1,strcat(measfilepath,measfilename,numstr,typestr,num2str(j+13),'.png'))
    
    intMeas = IntegrateWavelength(intsMeas, nuMeas, bands(j,1), bands(j,2));
    intModel = IntegrateWavelength(intsModel, nuModel, bands(j,1), bands(j,2));
    intConv = IntegrateWavelength(intsConv, nuModel, bands(j,1), bands(j,2));
    intModelH = IntegrateWavelength(intsModelH, nuModel, bands(j,1), bands(j,2));
    intConvH = IntegrateWavelength(intsConvH, nuModel, bands(j,1), bands(j,2));
    intModelC = IntegrateWavelength(intsModelC, nuModel, bands(j,1), bands(j,2));
    intConvC = IntegrateWavelength(intsConvC, nuModel, bands(j,1), bands(j,2));
  
    intsAll(j,:) = [intMeas, intModel, intConv, intModelH, intConvH, intModelC, intConvC];
    
end

end



