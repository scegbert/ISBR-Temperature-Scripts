function [Cv] = MvToCv(nu, mvs, coeffs, mus)

% Author: John Tobiasson
% Date: 07/28/2016
% Purpose: This function will take in a wavenumber (cm^-1) and an Mv and 
% convert the Mv to a Cv, with the intent of then converting the Mv to an 
% intensity with units of W/m^2/sr/um. This is done rather than converting
% directly to an intensity in this function in order to preserve the
% general process previously used.

% RevB note: This revision uses full band fits for A and B, but breaks down C
% into sections, so the overall equation becomes Cv = A = (x1*eta + x2) * 
% Mv^2 + (y1*eta + y2) * Mv + C, C being a function of wavenumber (eta) and
% it's order varying based on the band in which eta is located

% revC note: this revision splits the right section into two parts,
% necessary for the lower noise sets with more scans.

% RevD note: This revision is for use with different correlation values for
% each wavenumber. I'll find the closest wavenumber to what is in the list.

% Convert Mv to Cv
% find start and end rows

        
        for i = 1:1:size(nu, 1)                        
            % convert to Cvs
            Cv(i,:) = polyval(coeffs(i, 2:end), mvs(i), [], mus(i, 2:3));

        end
    
end