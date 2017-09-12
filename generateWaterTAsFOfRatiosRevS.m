%function [] = generateWaterTAsFOfRatiosRevS()

clear
clc
tic

% Creates a series of nested polynomial fit lines such that
% Tgasoptical = f(PL, f(P, f(Yh2o, f(E/A)))).
% Each f() sequence represents an additional layer of curve fits.
% In the lowest curve fit, the coefficients are stored. For all others,
% the resultant value, y = f(...), yields the coefficient for the next
% polynomial curve fit up in the ladder.

% note that under current conditions, code requires 19 hours to execute

% Cabs values have units of cm^2/molecule
% kappa values have units of 1/m

% Define band beginning and end points
% aBandStart = 5185;
% bBandStart = 5310;
% cBandStart = 5435;
% cBandEnd = 5560;
% eBandStart = 5615;
% eBandEnd = 5715;

aBandStart = 3800;
bBandStart = 3900;
cBandStart = 4000;
cBandEnd = 4100;
eBandStart = 3250;
eBandEnd = 3350;

% define constants
h1 = 6.62606957 * (10^(-34));
k1 = 1.3806488 * (10^(-23));
c0 = 2.99792458 * (10^8);

pathlength = .5; %(.02:.02:1)';
for pathlengthIter = 1:length(pathlength)
    
    pressure = 1; %[.5, 1, 2, 4, 8, 15, 30, 40]';
    for pressureIter = 1:length(pressure)
        
        yh2o = .2; %[.05, .10, .20, .30, .40]';
        for yh2oIter = 1:length(yh2o)
            
            temp = (300:100:3000);
            for tempIter = 1:length(temp)
                %print the conditions being calculated for reference
                PLpressureYH2O = [pathlength(pathlengthIter); pressure(pressureIter); yh2o(yh2oIter)]
                temperature = temp(tempIter)
                
                % calculate H2O kappas
                kappasH2O = calcGasKappasByEtaRevB('h2o', pressure(pressureIter),...
                    temp(tempIter), yh2o(yh2oIter), eBandStart, cBandEnd); %returns [nu, kappa]
                %aBandStart, eBandEnd); %returns [nu, kappa]
                
                % place holder to include CO2 kappas if desierd. Largely negligable for now
                % kappasH2O = calcGasKappasByEtaRevB('co2', pressure, temp, yco2, ...
                %    aBandStart, eBandEnd);
                
                % calculate emisivities
                emissivities = 1 - exp(-kappasH2O(:, 2) * pathlength(pathlengthIter));
                
                % calculate planck intensities
                planckInts = (2 * h1 * (c0^2) .* (kappasH2O(:, 1) .^5) * (10^4)) ...
                    ./ (exp((h1 * c0 .* kappasH2O(:, 1) * 100) / (k1 * temp(tempIter))) - 1);
                
                ints = emissivities .* planckInts;
                
                % find relevant wavenumbers in new sections
                [~, rowLowA] = min(abs(kappasH2O(:,1)-aBandStart));
                [~, rowLowB] = min(abs(kappasH2O(:,1)-bBandStart));
                [~, rowLowC] = min(abs(kappasH2O(:,1)-cBandStart));
                [~, rowHighC] = min(abs(kappasH2O(:,1)-cBandEnd));
                [~, rowLowE] = min(abs(kappasH2O(:,1)-eBandStart));
                [~, rowHighE] = min(abs(kappasH2O(:,1)-eBandEnd));
                
                % integrate bands. Flip both sets vertically to make the
                % integral positive
                aInt = trapz(flipud(10000./kappasH2O(rowLowA:rowLowB, 1)), flipud(ints(rowLowA:rowLowB)));
                bInt = trapz(flipud(10000./kappasH2O(rowLowB:rowLowC, 1)), flipud(ints(rowLowB:rowLowC)));
                cInt = trapz(flipud(10000./kappasH2O(rowLowC:rowHighC, 1)), flipud(ints(rowLowC:rowHighC)));
                eInt = trapz(flipud(10000./kappasH2O(rowLowE:rowHighE, 1)), flipud(ints(rowLowE:rowHighE)));
                
                % ratio intensities
                ratioEA(tempIter) = eInt / aInt;
                ratioEB(tempIter) = eInt / bInt;
                ratioEC(tempIter) = eInt / cInt;
                
                ratioCE(tempIter) = cInt / eInt;
                ratioCA(tempIter) = cInt / aInt;
                ratioBA(tempIter) = bInt / aInt;
                
                
%                 allEA(pathlengthIter, pressureIter, yh2oIter, tempIter) = ratioEA(tempIter);
%                 allEB(pathlengthIter, pressureIter, yh2oIter, tempIter) = ratioEB(tempIter);
%                 allEC(pathlengthIter, pressureIter, yh2oIter, tempIter) = ratioEC(tempIter);
                
                allEA(yh2oIter, tempIter) = ratioEA(tempIter);
                allEB(yh2oIter, tempIter) = ratioEB(tempIter);
                allEC(yh2oIter, tempIter) = ratioEC(tempIter);
                
                
                
            end %end temp loop
            
                            hold on
                plot(temp, ratioEA, 'LineWidth',2)
                plot(temp, ratioEB, 'LineWidth',2)
                plot(temp, ratioCE, 'LineWidth',2)
                plot(temp, ratioBA, 'LineWidth',2)
                plot(temp, ratioCA, 'LineWidth',2)
                legend('E/A','E/B','C/E','B/A','C/A')

            
            
            polyTempOrder = 3;
            
            polyTempA(yh2oIter,:) = polyfit(ratioEA, temp, polyTempOrder+2); %first polyfits
            polyTempB(yh2oIter,:) = polyfit(ratioEB, temp, polyTempOrder);
            polyTempC(yh2oIter,:) = polyfit(ratioEC, temp, polyTempOrder);
            
            tempCalcA = polyval(polyTempA(yh2oIter,:), ratioEA); 
            tempCalcB = polyval(polyTempB(yh2oIter,:), ratioEB);
            tempCalcC = polyval(polyTempC(yh2oIter,:), ratioEC);
     
%             plot(ratioEA, temp, ratioEB, temp, ratioEC, temp)
%             plot(temp, tempCalcA, temp, tempCalcB, temp, tempCalcC)
%             legend('A','B','C')
            
            clear ratioEA ratioEB ratioEC
            
        end %end yh2o loop
        
        polyYh2oOrder = 3;
        
        for i = 1:polyTempOrder+1 %once for each coefficient
            
            polyYh2oA(pressureIter,i,:) = polyfit(yh2o, polyTempA(:,i), polyYh2oOrder);
            polyYh2oB(pressureIter,i,:) = polyfit(yh2o, polyTempB(:,i), polyYh2oOrder);
            polyYh2oC(pressureIter,i,:) = polyfit(yh2o, polyTempC(:,i), polyYh2oOrder);
            
        end
        
        clear polyTempA polyTempB polyTempC
        
    end %end pressure loop
    
    polyPressureOrder = 3;
    
    
    for i = 1:polyTempOrder+1 %once for each coefficient
        for j = 1:polyYh2oOrder+1
            
            polyPressureA(pathlengthIter,i,j,:) = polyfit(pressure, polyYh2oA(:,i,j), polyPressureOrder);
            polyPressureB(pathlengthIter,i,j,:) = polyfit(pressure, polyYh2oB(:,i,j), polyPressureOrder);
            polyPressureC(pathlengthIter,i,j,:) = polyfit(pressure, polyYh2oC(:,i,j), polyPressureOrder);
            
        end
    end
    
    clear polyYh2oA polyYh2oB polyYh2oC
    
end %end pathlength loop

polyPathlengthOrder = 3;

for i = 1:polyTempOrder+1 %once for each coefficient
    for j = 1:polyYh2oOrder+1
        for k = 1:polyPressureOrder+1
            
            polyPathlengthA(i,j,k, :) = polyfit(pathlength, polyPressureA(:,i,j,k), polyPathlengthOrder);
            polyPathlengthB(i,j,k, :) = polyfit(pathlength, polyPressureB(:,i,j,k), polyPathlengthOrder);
            polyPathlengthC(i,j,k, :) = polyfit(pathlength, polyPressureC(:,i,j,k), polyPathlengthOrder);
            
        end
    end
end

clear polyPressureA polyPressureB polyPressureC

save('C:\Users\Scott\Documents\MATLAB\CabsDatabaseFiles\polyfit2', 'polyPathlengthA', 'polyPathlengthB', 'polyPathlengthC')
save('C:\Users\Scott\Documents\MATLAB\CabsDatabaseFiles\alldata', 'allEA', 'allEB', 'allEC')

toc


% organize polyPathlengths such that temperature can be calculated from them
% I think a similar algorithm will be needed to pull temperature out as was
% required to put it in there.
% would it be wise to look into using polyfitn at this point?
