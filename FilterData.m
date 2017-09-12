function [] = FilterData(calfilepath, calfilename, filterfilename, datafilepath,...
    starttemp, endtemp, tempstep, leftStart, rightEnd)

%% attempt at selective filtering based on absorption cross section. Work in progress

close all

%load Cabs coefficients, extract the data from a structure
H2Oraw = load(strcat(datafilepath,'relevantH2OCabsDatabaseFiles\Cabs h2o 300K Y1 P1'));
fieldName = fields(H2Oraw);
H2Odata = H2Oraw.(fieldName{1});

CO2raw = load(strcat(datafilepath,'relevantCO2CabsDatabaseFiles\Cabs co2 300K Y0 P1'));
fieldName = fields(CO2raw);
CO2data = CO2raw.(fieldName{1});

[~, leftBand] = min(abs(H2Odata(:,1) - leftStart));
[~, rightBand] = min(abs(H2Odata(:,1) - rightEnd));

%estimate H2O and CO2 concentration
ppmH2O = 7223.2;
ppmCO2 = 400; 

H2Oabsorb = H2Odata(leftBand:rightBand,2) * ppmH2O;
CO2absorb = CO2data(leftBand:rightBand,2) * ppmCO2;
nuAbsorb = H2Odata(leftBand:rightBand,1);

allAbsorb = H2Oabsorb + CO2absorb;

% plot(nuAbsorb,allAbsorb)

clear H2Oraw CO2raw H2Odata CO2data fieldname H2Oabsorb CO2absorb

%establish minimum value
[~, cutoffSplit1] = min(abs(nuAbsorb - 4500));
[~, cutoffSplit2] = min(abs(nuAbsorb - 6500));

cutoff(1:cutoffSplit1,1) = 1E-17; %units unknown, selected by looking at plot and filtered data
cutoff(cutoffSplit1:cutoffSplit2,1) = 2E-18; %best performance found by applying 2 different
cutoff(cutoffSplit2:length(allAbsorb),1) = 5E-18;

%create arrays for finding absorption rising and falling past cutoff
up = allAbsorb - cutoff;
up(up<0) = 0;
down = cutoff - allAbsorb;
down(down<0) = 0;

k = 1; %index in spectral array
j = 0; %index in absorption bands array
kFinal = length (up);

while k < kFinal && any(up) == 1
    j = j+1;
  
    absorbIndex(j,1) = find(up,1); %find section where cutoff starts
    k = absorbIndex(j,1); %shift spectral array index to new position
    
    down(1:k) = 0; %remove previously identified peaks
    
    absorbIndex(j,2) = find(down,1); %find section where cutoff ends
    k = absorbIndex(j,2); %shift spectral array index to new position
 
    up(1:k) = 0; %remove previously identified peaks

    absorbBands(j,:) = [nuAbsorb(absorbIndex(j,1)),nuAbsorb(absorbIndex(j,2))];
end

clear up down allAbsorb nuAbsorb

i = 0; %initialize counter variable
for T = starttemp:tempstep:endtemp
    tic
    i = i+1;    
    
%     if i == 2 || i == 5 || i == 10 || i == 13
%         i = i +1;
%     end               %for february 2016 calibration
    
    if i<10 %matches file format "...0001.csv" from OMNIC default
        numstr = strcat('000',num2str(i));
    elseif i<100
        numstr = strcat('00',num2str(i));
    elseif i<1000
        numstr = strcat('0',num2str(i));
    else
        numstr = num2str(i);
    end
    
    tempstr = num2str(T); %temperature value as a string

    % Load data
    data = csvread(strcat(calfilepath,calfilename,numstr,'.csv'));

    % Extract variables, prepare for filtering
    V = data(:,2); % Voltage
    nu = data(:,1); % Wave number
    clear data;
    
%% REMOVE SPECTRAL ABSORPTION

j = 0; %index of absorption band
Vabsorb = V; %copy over all data from original measurements
for j = 1:1:length(absorbBands)

    [~, startBand] = min(abs(nu-absorbBands(j,1)));
    [~, stopBand] = min(abs(nu-absorbBands(j,2)));
    
    %grab data outside of absorption range
    %average Vref to reduce impact of noise
    Vref = [mean(V(startBand - 5:startBand - 1)), mean(V(stopBand + 1:stopBand + 5))];
    nuRef = [nu(startBand - 1), nu(stopBand + 1)];
    
    Vabsorb(startBand:stopBand) = interp1(nuRef, Vref, nu(startBand:stopBand));
        
end
   

%% BUTTERWORTH SMOOTHING FILTER

    % Filter order
    nf = 5;

    % Filter cutoff frequency -- increasing this number decreases the
    % smoothing; decreasing this number increases the amount of smoothing
    fc = 0.004; %0.005; from Dr. Colton

    % Create filter
    [B,A] = butter(nf, fc);

    % Apply the butterworth filter to the previously filtered Vf data
    Vfilter = filtfilt(B, A, Vabsorb);

    % Plot unfiltered data
    figure(1);
    hold all
    plot(nu,V,'b');
    title(strcat(tempstr,'C Calibration Data'))
    xlabel('wavenumber(cm-1)')
    ylabel('Intensity (volts)')
    
    % Add "maximum" filtered data to plot
    plot(nu,Vabsorb,'g');
    
    % Add butterworth filtered data to plot
    figure(1);
    plot(nu,Vfilter,'r');
    legend('raw','maximum','smoothed')
    
    save(strcat(calfilepath,filterfilename,tempstr,'.mat'), 'nu',  'Vfilter');
 
    FilterDataTime = toc

end
end

%% ----------------Documentation-------------------
%{

For the 04282016 calibration

cutoff(1:cutoffSplit1,1) = 1E-17;
cutoff(cutoffSplit1:cutoffSplit2,1) = 2E-18; 
cutoff(cutoffSplit2:length(allAbsorb),1) = 5E-18;

fc = 0.004; 

For the 02172016 calibration

cutoff(1:cutoffSplit1,1) = 1E-17;
cutoff(cutoffSplit1:cutoffSplit2,1) = 1E-18; 
cutoff(cutoffSplit2:length(allAbsorb),1) = 5E-18;

fc = 0.005;



%}
