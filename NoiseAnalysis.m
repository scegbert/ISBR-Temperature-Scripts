function [overallSignalToNoise, waterToBB] = NoiseAnalysis(rawData, ...
            measGasInts, backgroundEmissivity, ...
            backgroundTemp, wavenumbers)

% check for needed wavenumbers for noise
if rawData(1, 1) <= 3600.1
    
    % can use 3600 - 3700
    [~, rowNoiseLow] = min(abs(rawData(:, 1)-3600));
    [~, rowNoiseHigh] = min(abs(rawData(:, 1)-3700));

else
    
    % must use 5800 - 5900
    [~, rowNoiseLow] = min(abs(rawData(:, 1)-5800));
    [~, rowNoiseHigh] = min(abs(rawData(:, 1)-5900));
    
end

% get signal rows - the signals seem to always have a BB peak in these
% areas

[~, rowSignalLow] = min(abs(rawData(:, 1)-4600));
[~, rowSignalHigh] = min(abs(rawData(:, 1)-4620));

% Get water rows
[~, rowWaterLow] = min(abs(wavenumbers(:, 1)-5350));
[~, rowWaterHigh] = min(abs(wavenumbers(:, 1)-5370));

% get max peak

maxSignalPeak = max(rawData(:,2));

% separate out arrays for pertinent data
noise = rawData(rowNoiseLow:rowNoiseHigh, 2);
% signal = rawData(rowSignalLow:rowSignalHigh, 2);
water = measGasInts(rowWaterLow:rowWaterHigh);

noiseHeight = mean(noise);

overallSignalToNoise = maxSignalPeak / noiseHeight;

% compare (peak) water intensity to peak BB intensity

    step = 0;
    for i=1:2:20
        step = step + 1/rawData(i, 1) - 1/rawData(i+1, 1);
    end
    avStep = (10^4)*step/10; 

    % Get water intensity
    
        waterSum = 0;
        
        k = 1;
        
        for j = rowWaterLow:1:rowWaterHigh 

            waterSum = water(k) + waterSum;
            
            k = k + 1;
            
        end
        
        avWaterInt = avStep .* waterSum ./ ( (10^4) ./ rawData(rowWaterLow, 1) - ...
                                            (10^4) ./ rawData(rowWaterHigh, 1) );
        
    % Get BB intensity
    
        bbSum = 0;
        
        h = 6.62606957 * (10^(-34));
        k = 1.3806488 * (10^(-23));
        c0 = 2.99792458 * (10^8);
        
        for m = rowWaterLow:1:rowWaterHigh 
                
            bbSum = bbSum + (2 * h * (c0^2) * ...
                (rawData(m, 1)^5) * (10^4)) / ( exp((h * c0 * ...
                rawData(m, 1) * 100) / (k * backgroundTemp)) - 1 ) * backgroundEmissivity;     
            
        end
        
        AvBBInt = (avStep .* bbSum) ./ ( (10^4) ./ rawData(rowWaterLow, 1) - ...
                                            (10^4) ./ rawData(rowWaterHigh, 1) );
        
    % Calculate waterToBB
    
    waterToBB = avWaterInt ./ AvBBInt;
    
end