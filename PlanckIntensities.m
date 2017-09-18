function [bbInts] = PlanckIntensities(temp, nu) %, waveType)

% Date: 07/15/2016
% Author: John Tobiasson
% Purpose: This function will calculate blackbody intensities for various
% temperatures for the wavenumbers of interest. The intensities will be as
% W/m^2/sr/um, based off of wavenumbers. This will be used as part of
% generating Cv values as functions of Mv and eta both.

% Temperatures must be in degrees C, and incrementable by 100 C.

% RevB Note: This revision no longer saves files, but is meant to be called

h = 6.62606957 * (10^(-34));
k = 1.3806488 * (10^(-23));
c0 = 2.99792458 * (10^8);

waveType = 'number';

if strcmp(waveType, 'number') == 1
    
    bbInts = (2 * h * (c0^2) .* (nu.^5) * (10^4)) ./ ...
        (exp((h * c0 .* nu * 100) / (k * temp)) - 1);
    
elseif strcmp(waveType, 'length') == 1
    
    bbInts = (2 * h * (c0^2) .* (nu.^3) * (10^8)) ./ ...
        (exp((h * c0 .* nu * 100) / (k * temp)) - 1);
    
else
    
    waveType = throwerror; %unkown waveType
    
end

end


