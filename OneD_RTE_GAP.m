function [intsModel, nuModel] = OneD_RTE_GAP( ...
                                    gasT, gasConcH2O, gasConcCO2, partT, partKappa,...
                                    wallT, wallEm, sensT, sensEm, ...
                                    pressure, pathLength, ...
                                    startWavenumber, stopWavenumber)
% Routine to calculate one-directional radiation intensity between
% two surfaces with a participating media.
% Based on an integration of the one-dimensional RTE.
% Based on Sections 8.2 - 8.7 of Modest "Radiative Heat Transfer," 1993.
%
% Intensity calculated at wall, passed cell-by-cell through medium.
% Intensity is attenuated and augmented based on uniform absorption and 
% scattering properties in each cell. 
% This formulation assumes constant, uniform cell properties, no
% in-scattering, only augmentation by cell emission. Also assumes only one
% temperature for cell, no differentiation between gas and particle temps.
%
% Written August 2016 by Brad Adams
%
% Modified Aug 2016 - Apr 2017 by John Tobiasson. This was modified to include
% comments explaining names of variables, and to modifiy this for use on a 
% spectral basis.
% Further modifications have expanded the ability of the code to include
% CO2 and particle radiative participation
%

%
% Constants and calculation tolerance
%
pi = 3.14159;
%sig = 5.67e-08; % Stefan-Boltzman constant
qtol = 0.0001;
    % Radiation constants
    h1 = 6.62606957 * (10^(-34));
    k1 = 1.3806488 * (10^(-23));
    c0 = 2.99792458 * (10^8);
    % wavenumbers of interest
nuModel = transpose(startWavenumber:0.005:stopWavenumber);
    
% Geometry Specification
% Assumes 1-D geom of length xtot divided into ni sections; 
% cell size (x-dimension) is specified by user.
% Current geometry is based on the BYU BFR: 0.75 m diameter, 2.4 m length.
%
% Sensor is assumed at i=ni or xf = 0.75 location.
%
xtot = pathLength;
ni = 10;
%
% Assume uniform spacing of dx for all cells now. 
%
xdx = xtot/ni;
for i=1:ni
    dx(i) = xdx;
end
%
% Surface and media properties - from NG-air testing in BFR; SI units
%
twni = sensT; % sensor temperature (K)
twi1 = wallT; % furnace wall temperature (K)
emwni = sensEm; % sensor emissivity
emwi1 = wallEm; % furnace wall emissivity
if size(gasT, 1) == ni % if one gas temp is input, the gas in all cells has
                       % the same temp. If an array of gas temps is passed,
                       % then each cell can have a different gas temp. CO2
                       % and H2O are given the same temp as each other.
    tg = gasT;
elseif size(gasT, 1) == 1
    tg = ones(ni, 1) .* gasT;
else
    tg = 0;
end
if size(partT, 1) == ni % particle temsp work as gas temps above.
    tp = partT;
elseif size(partT, 1) == 1
    tp = ones(ni, 1) .* partT;
else
    tp = 0;
end

% 
% Can add function calls to calculate gas absorption, soot absorption,
% particle absorption and particle scattering coefficients for each cell.
%
% This has been added for H2O and CO2
%
% gas absorption coefficient
    abskgH2O = zeros(size(nuModel, 1), ni); % define the H2O absorption
                                                % coefficient matrix. This
                                                % is organized as one
                                                % column for each cell,
                                                % with each row being for
                                                % each wavenumber.
    abskgCO2 = zeros(size(nuModel, 1), ni); % define the CO2 absorption
                                                % coefficient matrix
    [~, rowStart] = min(abs(nuModel-startWavenumber)); % find the start
                    % and stop row locations for the wavenumbers of
                    % interest. They will be the first and last rows here,
                    % but a modification of some sort could be made.
    [~, rowStop] = min(abs(nuModel-stopWavenumber));
    % load files with absorption coefficients
    for j = 1:1:ni
        
        if j > 1 % For subsequent iterations, only calculate new gas 
                 % absorption coefficients if the temperature are
                 % different. H2O and CO2 volumetric concentrations are
                 % assumed to be uniform throughout the path length for
                 % now.
            if tg(j) == tg(j - 1)
                abskgH2O(:, j) = abskgH2O(:, j - 1);
                abskgCO2(:, j) = abskgCO2(:, j - 1);
            else
               kappasH2O = calcGasKappasByEtaRevB('h2o', pressure, gasT(j), gasConcH2O, ...
                    startWavenumber, stopWavenumber);                                      
            
                kappasCO2 = calcGasKappasByEtaRevB('co2', pressure, gasT(j), gasConcCO2, ...
                        startWavenumber, 0);    

                abskgH2O(:, j) = kappasH2O(rowStart:rowStop, 2);

                try
                    abskgCO2(:, j) = kappasCO2(rowStart:rowStop, 2);
                catch
                    abskgCO2(1:size(kappasCO2,1), j) = kappasCO2(rowStart:end, 2);
%                     h2oSize = size(abskgH2O,1);
%                     co2Size = size(abskgCO2,1);
%                     zerosToAdd = zeros(h2oSize - co2Size,1);
%                     abskgCO2 = cat(1, abskgCO2, zerosToAdd);
                end
               
            end
            
        else
            
            kappasH2O = calcGasKappasByEtaRevB('h2o', pressure, gasT(j), gasConcH2O, ...
                    startWavenumber, stopWavenumber);                                      
            
            kappasCO2 = calcGasKappasByEtaRevB('co2', pressure, gasT(j), gasConcCO2, ...
                    startWavenumber, 0);
                
            abskgH2O(:, j) = kappasH2O(rowStart:rowStop, 2);
            
            try
                abskgCO2(:, j) = kappasCO2(rowStart:rowStop, 2);
            catch
                abskgCO2(1:size(kappasCO2,1), j) = kappasCO2(rowStart:end, 2);
            end
        end
    end
% Check the size of the CO2 array against that of the H2O array. Add zeros
% to fill out the CO2 array as necessary. This is necessary because
% Pearson's database does not provide CO2 absorption coefficients for the
% same number of wavenumbers as it does for H2O. CO2 arrays differ in size
% between sets of gas properties.
h2oSize = size(abskgH2O,1);
co2Size = size(abskgCO2,1);
if co2Size == h2oSize
else
    % assume co2 array is smaller than h2o (should always be the case, or
    % will throw an error later)
    zerosToAdd = zeros(h2oSize - co2Size,1);
    abskgCO2 = cat(1, abskgCO2, zerosToAdd);

end    
% absks = ones(size(wavenumbers, 1), ni).* 0; % soot absorption coefficient
abskp = ones(size(nuModel, 1), ni) .* partKappa; % particle absorption coefficient
scakp = ones(size(nuModel, 1), ni) .* 0; % particle scattering coefficient
%
% Calculate extinction coefficient, optical distance, and blackbody 
% emission intensity. 
%
extk = zeros(size(nuModel, 1), ni);
emtaup = zeros(size(nuModel, 1), ni);
emtaugH2O = zeros(size(nuModel, 1), ni);
emtaugCO2 = zeros(size(nuModel, 1), ni);
emtaus = zeros(size(nuModel, 1), ni);
abstau = zeros(size(nuModel, 1), ni);
bbintg = zeros(size(nuModel, 1), ni);
bbintp = zeros(size(nuModel, 1), ni);
for i=1:ni
    % extinction coefficient
%     extk(:, i) = abskgH2O(:, i) + abskgCO2(:, i) + absks(:, i) + abskp(:, i) + scakp(:, i);  
    extk(:, i) = abskgH2O(:, i) + abskgCO2(:, i) + abskp(:, i) + scakp(:, i); 
    % optical thickness                                              
    %tau(:, i) = extk(:, i) .* dx(i); % this is not currently used in the
    % code.
    % emissive particle optical thickness
    emtaup(:, i) = abskp(:, i) .* dx(i);
    % emissive gas optical thickness
    emtaugH2O(:, i) = abskgH2O(:, i) .* dx(i);
    emtaugCO2(:, i) = abskgCO2(:, i) .* dx(i);
    % emissive soot optical thickness
%     emtaus(:, i) = absks(:, i) .* dx(i);
    % absorptive optical thickness
    abstau(:, i) = extk(:, i) .* dx(i);
    % blackbody intensity
    bbintg(:, i) = (2 * h1 * (c0^2) .* (nuModel.^5) * ...
                (10^4)) ./ ( exp((h1 * c0 .* nuModel * 100) / ...
                (k1 * tg(i))) - 1 ); % Calcualted for each cell, based on
                                     % the gas temperature
    bbintp(:, i) = (2 * h1 * (c0^2) .* (nuModel.^5) * ...
                (10^4)) ./ ( exp((h1 * c0 .* nuModel * 100) / ...
                (k1 * tp(i))) - 1 ); % particle black body intensity for
                                     % each cell           
end
%
% Initialize boundary incident fluxes.
%
qwni = zeros(size(nuModel, 1), 1); % incident flux at sensor
%qwi1 = zeros(size(wavenumbers, 1), 1); % incident flux furnace wall
qwnst = ones(size(nuModel, 1), 1) .* 1e-6; % stored flux at the wall,
            % used to check for convergence.
qw1st = ones(size(nuModel, 1), 1) .* 1e-6; % stored flus at the sensor,
            % also used to check for convergence.
xintbc = zeros(size(nuModel, 1), 10); % incident flux for each cell
% 
% Start sweep from sensor side, i=ni.
%
while 1
  for i=ni:-1:1
    if (i == ni)        
        % intensity (in x direction) at i cell face
            xintbc(:, i) = emwni * ((2 * h1 * (c0^2) .* (nuModel.^5) * ...
                (10^4)) ./ ( exp((h1 * c0 .* nuModel * 100) / ...
                (k1 * twni)) - 1 )) + (1-emwni) .* qwni ./ pi;
    else
        xintbc(:, i) = xintfc(:, i+1);
    end
    % intensity (in x direction) at i + 1 cell face
    xintfc(:, i) = xintbc(:, i).*exp(-abstau(:, i)) ...
                   + (1-exp(-(abskp(:,i)+abskgH2O(:,i)+abskgCO2(:,i)).*dx(i))) ...
                   ./ (abskp(:, i) + abskgH2O(:, i) + abskgCO2(:,i)) ...
                   .*(bbintg(:, i).* abskgH2O(:,i) + bbintg(:, i).* abskgCO2(:,i)...
                   + bbintp(:, i).*abskp(:,i));
  end
  qwi1 = xintfc(:, 1) * pi;
%
% Sweep from furnace wall side, i=1
%
  for i=1:ni
    if (i == 1)
            xintbc(:,i) = emwi1 * ((2 * h1 * (c0^2) .* (nuModel.^5) * ...
                (10^4)) ./ ( exp((h1 * c0 .* nuModel * 100) / ...
                (k1 * twi1)) - 1 )) + (emwi1) .* qwi1 ./ pi;
            % This is based on the idea that a cavity, not a wall, is emitting and reflecting here.
            % (1 - emwi1) .* qwi1 ./ pi; for a wall that is emitting and reflecting
            % (emwi1) .* qwi1 ./ pi; for a cavity
    else
        xintbc(:,i) = xintfc(:,i-1);
    end
    xintfc(:, i) = xintbc(:, i).*exp(-abstau(:, i)) ...
                   + (1-exp(-(abskp(:,i)+abskgH2O(:,i)+abskgCO2(:,i)).*dx(i))) ...
                   ./ (abskp(:, i) + abskgH2O(:, i) + abskgCO2(:,i)) ...
                   .*(bbintg(:, i).* abskgH2O(:,i) + bbintg(:, i).* abskgCO2(:,i)...
                   + bbintp(:, i).*abskp(:,i));
  end
  qwni = xintfc(:,ni) * pi;
  wallint = xintbc(:,1);
  sensint = xintfc(:,ni);
intsModel = xintfc(:, ni);
%
% Check for convergence.
%
  qnerr = (qwni - qwnst) ./ qwnst;
  q1err = (qwi1 - qw1st) ./ qw1st;
  qwnst = qwni;
  qw1st = qwi1;
  qnerrA = mean(qnerr);
  q1errA = mean(q1err);
  if (qnerrA < qtol && q1errA < qtol)
     break
  end
end

end
