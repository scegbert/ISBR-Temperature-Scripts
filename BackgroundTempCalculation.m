function [tempBroadband, emBroadband, intLeft, intRight]...
    = BackgroundTempCalculation(bandsBroadband, nu, ints)

% This function will calculate the background temperature and
% emissivity based on 2-color pyrometry. This is specifically for use when
% processing data without particles.

%% filter measured data

%[intsMeasFiltered] = FilterBackground(datafilepath, BgBands, wavenumbers, measGasInts);

%% calculate background temperature
% Get row information
rowLowLeft = bandsBroadband(1);
rowHighLeft = bandsBroadband(2);
rowLowRight = bandsBroadband(3);
rowHighRight = bandsBroadband(4);

% Define constants (from online, similar but with more decimals ...
% than in my heat transfer book)
h = 6.62606957 * (10^(-34));
k = 1.3806488 * (10^(-23));
c0 = 2.99792458 * (10^8);

% Calculate the background intensities based on the background
% temperature guesses. Calcs Planck intensities and converts measured
% values

nuLeft = 10^4 ./ nu(rowLowLeft:rowHighLeft);
nuRight = 10^4 ./ nu(rowLowRight:rowHighRight);

intsLeft = ints(rowLowLeft:rowHighLeft);
intsRight = ints(rowLowRight:rowHighRight);

% Change sums to integrals
intLeft = trapz(flipud(nuLeft), flipud(intsLeft));
intRight = trapz(flipud(nuRight), flipud(intsRight));

if intLeft < .5 % if the signal is so weak we would get a bad value if we tried
    
    tempBroadband = 900;
    emBroadband = .0125;
    
else
    % ratio bands
    ratioMeasLR = intLeft / intRight;
    
    % calculate temperature by comparing to Planck curve (two color pyrometry method)
    tempRange = 200:100:3000; %less refined initial run
    
    for i = 1:2 %execute twice, once less refined and once more refined
        
        for j = 1:1:length(tempRange)
            
            intsPlanckLeft = (2 * h * (c0^2) * ((nu(rowLowLeft:rowHighLeft)).^5) .* ...
                (10^4)) ./ (exp((h .* c0 .* (nu(rowLowLeft:rowHighLeft)) .* 100) / ...
                (k .* tempRange(j))) - 1 );
            
            intsPlanckRight = (2 * h * (c0^2) * ((nu(rowLowRight:rowHighRight)).^5) .* ...
                (10^4)) ./ (exp((h .* c0 .* (nu(rowLowRight:rowHighRight)) .* 100) / ...
                (k .* tempRange(j))) - 1 );
            
            intPlanckLeft = trapz(flipud(nuLeft), flipud(intsPlanckLeft));
            intPlanckRight = trapz(flipud(nuRight), flipud(intsPlanckRight));
            
            RatioLRPlanck(j,1) = intPlanckLeft / intPlanckRight;
            
        end
        
        tempBroadband = interp1(RatioLRPlanck, tempRange, ratioMeasLR);
        
        tempRange = (tempBroadband-150):.5:(tempBroadband+150); %more refined second run
        
    end
    
    %% Get the emissivity of the background
    
    % Add up total planck values for em = 1, for each section
    intsPlanckLeft = (2 .* h .* (c0^2) .* ((nu(rowLowLeft:rowHighLeft)).^5) .* ...
        (10^4)) ./ ( exp((h .* c0 .* (nu(rowLowLeft:rowHighLeft)) .* 100) / ...
        (k .* tempBroadband)) - 1 );
    
    intsPlanckRight = (2 .* h .* (c0^2) .* ((nu(rowLowRight:rowHighRight)).^5) .* ...
        (10^4)) ./ ( exp((h .* c0 .* (nu(rowLowRight:rowHighRight)) .* 100) / ...
        (k .* tempBroadband)) - 1 );
    
    % convert planck sums to integrals for total intensity
    intPlanckLeft = trapz(flipud(nuLeft), flipud(intsPlanckLeft));
    intPlanckRight = trapz(flipud(nuRight), flipud(intsPlanckRight));
    
    emBroadbandLeft = intLeft / intPlanckLeft;
    emBroadbandRight = intRight / intPlanckRight;
     
    % get emissivity from intensities
    emBroadband = mean([emBroadbandLeft, emBroadbandRight]);
    
end

end