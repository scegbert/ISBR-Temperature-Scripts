function [bandint] = OneD_RTE_Gas( gasT, gasConcH2O, gasConcCO2, ...
                                   wallT, wallEm, sensT, sensEm, ...
                                   pressure, pathLength, ...
                                   startWavenumber, stopWavenumber)
%function [] = OneD_RTE_Gas
%--------------------------------------------------------------------------
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
% Modified Aug 2016 - Apr 2017 by John Tobiasson. This was modified to 
% include comments explaining names of variables, and to modifiy this for 
% use on a spectral basis.
% Further modifications have expanded the ability of the code to include
% CO2 and particle radiative participation.
%
% Modified May 2017 by Brad Adams to perform calculation of total intensity
% over all specified wavelengths.
%--------------------------------------------------------------------------
% Inputs passed in.
%
%gasT = 1387.;
%gasConcH2O = 0.168; 
%gasConcCO2 = 0.084;
%wallT = 300.;
%wallEm = 1.0; 
%sensT = 300.;
%sensEm = 1.0;
%pressure = 0.9;
%pathLength = 0.65;
%startWavenumber = 5000.;
%stopWavenumber = 5500.;
%
% Constants and calculation tolerance
%
pi = 3.14159;
sig = 5.67e-08; % Stefan-Boltzman constant
qtol = 0.0001;
%
% Radiation constants - with these constants must adjust blackbody 
% intensity calculation to produce desired units of W/m^2/um for wavelength
% (use conversion of 10^4 in numerator) or W/m^2/cm for wavenumber (use
% conversion of 10^8 in numerator)
%
h1 = 6.62607 * (10^(-34));
k1 = 1.38065 * (10^(-23));
c0 = 2.99792 * (10^8);
%
% wavenumbers of interest - arrange so wavenumbers are in column 1 of array
wavenumbers = transpose(startWavenumber:0.005:stopWavenumber);
intensityEta = 0.0;
%
% Geometry Specification
% Assumes 1-D geom of length xtot divided into ni sections; 
% cell size (x-dimension) is specified by user.
% 
% Code originally used BYU BFR geometry of 0.65 m diameter for pathlength.
% FTIR sensor taking spectral data in narrow waveband was assumed to be
% located at i=ni or xf = 0.65 location. Now read in pathlength and divide
% into 10 cells of similar composition and temperature.
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
emwi1 = wallEm; % furnace wall emissivity; may be cold cavity
partT = gasT; % set particle temps equal to gas temps
partKappa = 0.; % for now, set particle absorption to zero
%
% Fill gas and particle temperature arrays. 
%
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
if size(partT, 1) == ni % particle temps work as gas temps above.
    tp = partT;
elseif size(partT, 1) == 1
    tp = ones(ni, 1) .* partT;
else
    tp = 0;
end
% 
% Can add function calls to calculate gas absorption, soot absorption,
% particle absorption and particle scattering coefficients for each cell.
% For now, since in-scattering is not accounted for, do not include
% out-scattering either.
%
% This has been added for gas absorption from H2O and CO2
%
% gas absorption coefficient
%
abskgH2O = zeros(size(wavenumbers, 1), ni); % define the H2O absorption
                                            % coefficient matrix. This
                                            % is organized as one
                                            % column for each cell,
                                            % with each row being for
                                            % each wavenumber.
abskgCO2 = zeros(size(wavenumbers, 1), ni); % define the CO2 absorption
                                            % coefficient matrix
[~, rowStart] = min(abs(wavenumbers-startWavenumber)); % find the start
                    % and stop row locations for the wavenumbers of
                    % interest. They will be the first and last rows here,
                    % but a modification of some sort could be made.
[~, rowStop] = min(abs(wavenumbers-stopWavenumber));
%
% Calculate absorption coefficients based on spectral gas absorption data.
% Gas absorption data comes from routine calcGasKappasByEtaRevB for H2O
% and CO2 data. This routine in turn uses data developed by Pearson et al
% (J. Quant Spect Rad Trans 143 (2014) 100-110) and Pearson thesis based
% on the HITEMP 2010 spectroscopic database.
% The data files for H2O and CO2 are required for this code to work.
%
% Fill absorption coefficient arrays with calculated values.
%
for j = 1:1:ni
        
    if j > 1 % For subsequent iterations, only calculate new gas 
             % absorption coefficients if the temperatures are
             % different. H2O and CO2 volumetric concentrations are
             % assumed to be uniform throughout the path length for
             % now.
        if tg(j) == tg(j - 1)
            abskgH2O(:, j) = abskgH2O(:, j - 1);
            abskgCO2(:, j) = abskgCO2(:, j - 1);
        else
            kappasH2O = calcGasKappasByEtaRevB('h2o',pressure,gasT(j), ...
                        gasConcH2O,startWavenumber,stopWavenumber);                                      
            
            kappasCO2 = calcGasKappasByEtaRevB('co2',pressure,gasT(j), ...
                        gasConcCO2,startWavenumber, 0);    

            abskgH2O(:, j) = kappasH2O(rowStart:rowStop, 2);
            % Some CO2 arrays have more entries than others. Try
            % assuming they are the same size, if not, just add the 
            % full new array leaving remainder as zeroes
            try
              abskgCO2(:, j) = kappasCO2(rowStart:rowStop, 2);
            catch
              abskgCO2(1:size(kappasCO2,1), j) = kappasCO2(rowStart:end,2);
%                 h2oSize = size(abskgH2O,1);
%                 co2Size = size(abskgCO2,1);
%                 zerosToAdd = zeros(h2oSize - co2Size,1);
%                 abskgCO2 = cat(1, abskgCO2, zerosToAdd);
            end             
        end          
    else
%
% First iteration to calculate gas absorption coefficients (j=1)
%
        kappasH2O = calcGasKappasByEtaRevB('h2o', pressure, gasT(j), ...
                    gasConcH2O, startWavenumber, stopWavenumber);                                      
   
        kappasCO2 = calcGasKappasByEtaRevB('co2', pressure, gasT(j), ...
                    gasConcCO2, startWavenumber, 0);

        abskgH2O(:, j) = kappasH2O(rowStart:rowStop, 2);
            
        try
            abskgCO2(:, j) = kappasCO2(rowStart:rowStop, 2);
        catch
            abskgCO2(1:size(kappasCO2,1), j) = kappasCO2(rowStart:end, 2);
        end
    end
end
%
% Check the size of the CO2 array against that of the H2O array. Add zeros
% to fill out the CO2 array as necessary. This is necessary because
% Pearson's database does not provide CO2 absorption coefficients for the
% same number of wavenumbers as it does for H2O. CO2 arrays differ in size
% between sets of gas properties.
%
h2oSize = size(abskgH2O,1);
co2Size = size(abskgCO2,1);
if co2Size ~= h2oSize
    % assume co2 array is smaller than h2o (should always be the case, or
    % will throw an error later)
    zerosToAdd = zeros(h2oSize - co2Size,1);
    abskgCO2 = cat(1, abskgCO2, zerosToAdd);
end   
%
% Initialize soot and particle absorption and particle scattering
% coefficients. Note these are spectral arrays similar to kappasH2O.
% With no additional calculations for these properties at this time,
% these arrays will be filled with zeroes, except the partKappa value
% which could be specified by the user (others could also be specifed).
% In particular, since we do not consider in-scattering, do not consider
% out-scattering at this time.
%
absks = ones(size(wavenumbers, 1), ni).* 0; % soot absorption coefficient
abskp = ones(size(wavenumbers, 1), ni) .* partKappa; % particle absorption coefficient
%scakp = ones(size(wavenumbers, 1), ni) .* 0; % particle scattering coefficient
%
% Initialize then calculate extinction coefficient, optical thickness
% (defined as absorption or extinction coefficients times optical path
% length, which here is the cell width), and blackbody emission intensity
% arrays for all cells along pathlength. 
%
abskgt = zeros(size(wavenumbers, 1), ni);
extkg = zeros(size(wavenumbers, 1), ni);
tau = zeros(size(wavenumbers, 1), ni);
bbintg = zeros(size(wavenumbers, 1), ni);
bbintp = zeros(size(wavenumbers, 1), ni);
bintav = zeros(size(wavenumbers, 1), ni);
%
for i=1:ni
    % gas-related absorption coefficients (H2O, CO2, soot)
    abskgt(:, i) = abskgH2O(:, i) + abskgCO2(:, i) + absks(:, i);
    
    % extinction coefficient (gas, soot, particle absorption)
    % Because we do not consider in-scattering here, do not include 
    % particle out-scattering.
    % extkg(:, i) = abskgt(:, i) + abskp(:, i) + scakp(:, i);  
    extkg(:, i) = abskgt(:, i) + abskp(:, i);  
    
    % optical thickness                                              
    tau(:, i) = extkg(:, i) .* dx(i); 
    
    % Blackbody emission intensity - conversion factor in numerator depends
    % on which units for bb intensity. To produce desired units of W/m^2/um 
    % for wavelength use conversion of 10^24 in numerator, or for 
    % W/m^2/cm-1 for wavenumber use conversion of 10^8 in numerator.
    % Also need adjustment in exponent in denominator.
    % For now we will calculate all intensities in cm-1 (10^8 adjustment).
    % Include this in xintbc calcs as well.
    %
    % Gas and particle bb emission intensities calculated for each cell
    % based on local gas and particle temperature, respectively.
    % Also calculate an average bb intensity which includes gas and 
    % particle effects, for use in cell intensity calcs.
    %
%    bbi5185eta = (2 * h1 * (c0^2) * (5185.2^3) * (10^8)) ...
%               / ( exp((h1 * c0 * 5185.2 * 100) / (k1 * tg(i))) - 1 );
%    bbi5185mu = (2 * h1 * (c0^2) * (10^24)) / (1.9286^5) ...
%               / ( exp((h1 * c0 * 10^6) / (k1 * 1.9286 * tg(i))) - 1 );
    bbintg(:, i) = (2 * h1 * (c0^2) .* (wavenumbers.^3) * (10^8)) ...
                 ./ ( exp((h1 * c0 .* wavenumbers * 100) / ...
                 (k1 * tg(i))) - 1 ); 
      %
    bbintp(:, i) = (2 * h1 * (c0^2) .* (wavenumbers.^3) * (10^8)) ...
                 ./ ( exp((h1 * c0 .* wavenumbers * 100) / ...
                 (k1 * tp(i))) - 1 );
    %
    bintav(:, i) = (bbintg(:, i) .* abskgt(:, i) + ...
                    bbintp(:, i) .* abskp(:,i)) ...
                ./ (abskgt(:, i) + abskp(:, i));
end
%
% Initialize boundary incident fluxes.
%
qwni = zeros(size(wavenumbers, 1), 1); % incident flux at sensor
qwi1 = zeros(size(wavenumbers, 1), 1); % incident flux furnace wall/cavity
qwnst = ones(size(wavenumbers, 1), 1) .* 1e-6; % stored flux at the wall,
            % used to check for convergence.
qw1st = ones(size(wavenumbers, 1), 1) .* 1e-6; % stored flux at the sensor,
            % also used to check for convergence.
xintbc = zeros(size(wavenumbers, 1), 10); % incident flux for each cell
% 
% Assume furnace wall/cold cavity at i=1 and sensor at i=ni
%
% Start sweep from sensor side, i=ni.
%
while 1
  for i=ni:-1:1
    if (i == ni)        
        % Intensity in -x direction originating at i cell face
        % Intensity has units of W/m^2/sr/cm-1 (use 10^8 conversion factor)
        % Assume gray, diffuse surfaces so absorptivity = emissivity and
        % reflectivity = 1 - absorptivity,
        % refected intensity is then qincident*reflectivity/pi
        xintbc(:, i) = emwni * ((2 * h1 * (c0^2) .* (wavenumbers.^3) * ...
                (10^8)) ./ ( exp((h1 * c0 .* wavenumbers * 100) / ...
                (k1 * twni)) - 1 )) + (1-emwni) .* qwni ./ pi;
    else
        xintbc(:, i) = xintfc(:, i+1);
    end
    % Intensity (in -x direction) at i-1 cell face
    % Calculate change in intensity moving across cell i
    % xintfc has units of W/m^2/sr/cm (based on 10^8 conversion factor)
    %
    xintfc(:, i) = xintbc(:, i) .* exp(-tau(:, i)) ...
                 + bintav(:, i) .* (1 - exp(-tau(:, i)));
  end
  % Convert to incident heat flux assuming diffuse surface
  qwi1 = xintfc(:, 1) .* pi;
%
% Sweep from furnace wall/cold cavity side, i=1
%
  for i=1:ni
    if (i == 1)
        % Intensity in x direction originating at i cell face
        % Intensity has units of W/m^2/sr/cm-1 (use 10^8 conversion factor)
        % Assume gray, diffuse surfaces so absorptivity = emissivity and
        % reflectivity = 1 - absorptivity,
        % refected intensity is then qincident*reflectivity/pi
        % Here the absorptivity in the "(1-emwi1) .* qwi1 ./ pi" term 
        % has been replaced by a fixed value of 1-0.9583, which has been
        % tuned to match data from reactor based on the idea that a 
        % "composite" wall comprised of a mostly blackbody with a small
        % impact of surrounding physical wall is reflecting here.
        % This hard code makes the emissivity, emwi1, and reflectivity, 
        % 1-emwi1, the same for the "composite" wall. 
%        xintbc(:,i) = emwi1 * ((2 * h1 * (c0^2) .* (wavenumbers.^5) * ...
%            (10^8)) ./ ( exp((h1 * c0 .* wavenumbers * 100) / ...
%            (k1 * twi1)) - 1 )) + (1 - emwi1) .* qwi1 ./ pi;
        xintbc(:,i) = emwi1 * ((2 * h1 * (c0^2) .* (wavenumbers.^3) * ...
            (10^8)) ./ ( exp((h1 * c0 .* wavenumbers * 100) / ...
            (k1 * twi1)) - 1 )) + (emwi1) .* qwi1 ./ pi;
    else
        xintbc(:,i) = xintfc(:,i-1);
    end
    % Intensity (in x direction) at i+1 cell face
    % Calculate change in intensity moving across cell i
    % xintfc has units of W/m^2/sr/cm (based on 10^8 conversion factor)
    %
    xintfc(:, i) = xintbc(:, i) .* exp(-tau(:, i)) ...
                 + bintav(:, i) .* (1 - exp(-tau(:, i)));
  end
  % Convert to incident heat flux assuming diffuse surface
  qwni = xintfc(:,ni) .* pi; 
%
% intensityEta returns spectral array of intensities arriving at sensor
% for given wavenumber band
%
  intensityEta = xintfc(:, ni);
%
% Check for convergence.
%
  qnerr = (qwni - qwnst) ./ qwnst;
  q1err = (qwi1 - qw1st) ./ qw1st;
  qwnst = qwni;
  qw1st = qwi1;
  qnerrA = max(qnerr);
  q1errA = max(q1err);
  if (qnerrA < qtol && q1errA < qtol)
  %
  % Now calculate the total intensity for the specified wavenumber band.
  % Idea is to divide the total spectrum of interest (200 - 20,000 1/cm)
  % into 20 bands of 1000 1/cm each at increments of 0.005 1/cm.
  % This will avoid overrunning array sizes and give some prelim output.
  %
  % Assume intensity is constant across 0.005 1/cm bandwidth. Integrate
  % simply by summing over the bands in increments of 0.005 1/cm.
  % 
    bandint = 0.0;
    nbnd = (stopWavenumber - startWavenumber) / 0.005;
    for i=1:nbnd+1
      bandint = bandint + intensityEta(i) * 0.005;
    end
     break
  end
end

%--------------------------------------------------------------------------
% Plotting could happen, but has been disabled for now. It may need to be
% fixed up some to make it actually work, as I haven't used it for some
% time. - JRT 22 April 2017
% Plot data
%
%fint(1) = emwi1 * sig * twi1^4 / pi + (1-emwi1) * qwi1 / pi;
%wint(1) = emwi1 * sig * twi1^4 / pi + (1-emwi1) * qwi1 / pi;
%xf(1) = 0.0;
% for (i=2:ni+1)
%    wint(i) = wint(1);
%    fint(i) = xintfc(i-1);
%    xf(i) = xf(i-1) + dx(i-1);
% end
%
% figure (1);
% plot(xf,wint,xf,fint)
% xlabel('Distance (m)');
% ylabel('Intensity (W/m^2/sr)');
% legend('Wall Intensity','Gas Intensity','Location','northwest');
%ax = gca;
%ax.YAxis.TickLabelFormat = '%,8.1f';

% subplot(2,1,1);
% plot(xf,fint);
% xlabel('Distance (m)');
% ylabel('Intensity (W/m^2/sr)');
% legend('Gas Intensity');
% subplot(2,1,2);
% plot(xf,wint);
% xlabel('Distance (m)');
% ylabel('Intensity (W/m^2/sr)');
% legend('Wall Intensity');
%--------------------------------------------------------------------------
end
