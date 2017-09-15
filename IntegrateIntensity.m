function [int] = IntegrateWavelength(ints, nu, start, stop) %, waveType)

% integrates a function within the desired band
% integration units are with respect to wavelength

[~, indexStart] = min(abs(nu-start));
[~, indexStop] = min(abs(nu-stop));

% Get Cvs (correction constants), calc intensities, and sum bands
intsTruncated = ints(indexStart:indexStop);

waveType = 'length';

if indexStart == indexStop
    
    int = 0;
    
elseif strcmp(waveType, 'number') == 1
    
    % Change sums to integrals, flipping arrays so the integrals are +
    int = trapz(nu, intsTruncated);
    
elseif strcmp(waveType, 'length') == 1
    
    % get wavelength band values
    lambda = 10^4 ./ nu(indexStart:indexStop);
    
    % Change sums to integrals, flipping arrays so the integrals are +
    int = trapz(flipud(lambda), flipud(intsTruncated));
    
else
    
    waveType = throwerror; %unkown waveType
    
end

end
