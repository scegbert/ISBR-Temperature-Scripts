function [] = PlotCalibratedData(calfilepath, filterfilename, starttemp, endtemp, tempstep, leftStart, rightEnd)

load(strcat(calfilepath,'mvToCvCoeffs.mat'));

    for T = starttemp:tempstep:endtemp

    rawData = xlsread(strcat(calfilepath, filterfilename,num2str(T),'.xlsx'));
            
    rawWavenumbers = rawData(:,1);
    [~, rowLowLeft] = min(abs(rawWavenumbers-leftStart));
    [~, rowHighRight] = min(abs(rawWavenumbers-rightEnd));

    newWavenumbers = rawData(rowLowLeft:rowHighRight, 1);
    voltages = rawData(rowLowLeft:rowHighRight, 2);
    
    Cvs = MvToCv(newWavenumbers, voltages, coeffs, mus);
    measGasInts = (voltages ./ Cvs);    
    plot(newWavenumbers,measGasInts)
    ylabel('emission intensity')
    xlabel('wavenumber [cm-1]')
    
   pause
    
   end
   end