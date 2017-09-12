function [] = CalibrationCurves(calfilepath, filterfilename, tempStart,...
    tempEnd, tempStep, leftStart, rightEnd)

close all
tic

% set order of polynomials for mvs
mvsOrder = 4; %John used 3

% set low and high temperatures in degrees C
tempLow = min(tempStart,tempEnd);
tempHigh = max(tempStart,tempEnd);
tempSpacing = abs(tempStep);

tempstr = num2str(tempLow);
load(strcat(calfilepath,filterfilename,tempstr,'.mat'));

[~, rowLowLeft] = min(abs(nu-leftStart)); %index of start and stop wavenumbers
[~, rowHighRight] = min(abs(nu-rightEnd));
nuTruncated = nu(rowLowLeft:rowHighRight);

% read in original wavenumbers
i = 0;
clear nu Vfilter

for temp = tempLow:tempSpacing:tempHigh
    
    i = i+1;
    
    tempstr = num2str(temp);
    load(strcat(calfilepath,filterfilename,tempstr,'.mat'));
    % variables: 'nu',  'Vfilter'
    mvs(:,i) = Vf2(rowLowLeft:rowHighRight); %only desired wavenumbers are saved
    
    bbInts(:,i) = PlanckIntensities(temp, nuTruncated);
    
    clear nu Vfilter
    
end

% create cvs
cvs = mvs ./ bbInts;

% create Cv = f(Mv) for each eta
coeffs = zeros(size(nuTruncated, 1), (mvsOrder + 2));
mus = zeros(size(coeffs, 1), 3);
coeffs(:, 1) = nuTruncated;
mus(:, 1) = nuTruncated;

for m = 1:size(mus, 1)
    [fit, ~, mu] = polyfit(mvs(m, :), cvs(m, :), mvsOrder);
    coeffs(m, 2:end) = fit;
    mus(m, 2:3) = mu;
    
end

% plot coeffs to verify that they follow bb curve
i = 1;
hold all
for temp = tempLow:tempSpacing:tempHigh
    plot(nuTruncated, bbInts(:, i), 'o')
    title('Calibration Curve Fit')
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Intensity (W/{m^2}/sr/cm-1');
    
    Cvs(:,i) = MvToCv(nuTruncated, mvs(:,i), coeffs, mus);
    measInts(:,i) = (mvs(:,i) ./ Cvs(:,i));
    
    plot(nuTruncated, measInts(:,i),'r')
    legend('BB Intensity','Curve Fit')
    
    i = i+1;
end

% save coeffs
save(strcat(calfilepath,'mvToCvCoeffs.mat'), 'coeffs',  'mus');

% save groups of data
% xlswrite(strcat(calfilepath, 'wavenumbers.xlsx'), originalWavenumbers);
% xlswrite(strcat(calfilepath, 'mvs.xlsx'), mvs);
% xlswrite(strcat(calfilepath, 'cvs.xlsx'), cvs);
% xlswrite(strcat(calfilepath, 'bbInts.xlsx'), bbInts);

CalibrationCurveTime = toc
pause
end
