function [int] = IntegrateWavelength(ints, nu, start, stop)

% integrates a function within the desired band
% integration units are with respect to wavelength

[~, indexStart] = min(abs(nu-start));
[~, indexStop] = min(abs(nu-stop));

if indexStart == indexStop
    
    int = 0;
    
else
    
    % Get Cvs (correction constants), calc intensities, and sum bands
    intsTruncated = ints(indexStart:indexStop);
    
    % get wavelength band values
    lambda = 10^4 ./ nu(indexStart:indexStop);
    
    % Change sums to integrals, flipping arrays so the integrals are +
    int = trapz(flipud(lambda), flipud(intsTruncated));
end
end
