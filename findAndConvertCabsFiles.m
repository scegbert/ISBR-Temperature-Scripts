% This code will find relevant Cabs files in the database for water. Once a
% relevant file has been found, it will be loaded, split into wavenumbers
% and Cabs values, and saved as a .mat file. These new files will load in
% .5s, as compared to 6.5-25s with various methods of loading the currently
% existing files.
%
% Author: John Tobiasson
% Date: 08/09/2016

% clear;
% clc;

% search through the relevant folder
files = dir('D:\Cabs_Absorption Coefficients_h2o\*.txt');    
for file = files'
    % make sure it's a .txt file
    if strcmp(file.name(end-2:end), 'txt') == 1
        % make sure it's a file for the right gas
        if strcmp(file.name(6:8), 'h2o') == 1
            % search through it looking for the temperature
                % find K
                tempLoc = strfind(file.name, 'K');
            tempStr = file.name(10:tempLoc-1);
            tempNum = str2double(tempStr);

            if tempNum == 500 %>= 299 && tempNum <= 700            
            % if it's a relevant temperature, find the concentration
                concStart = tempLoc + 3;
                postConcBlank = strfind(file.name(concStart:end), ' ');
                concEnd = (postConcBlank - 2) + concStart;
                concStr = file.name(concStart:concEnd);

                if strcmp(concStr, '0') == 1 || strcmp(concStr, '05') == 1 || strcmp(concStr, '1') == 1  ...
                     || strcmp(concStr, '2') == 1 || strcmp(concStr, '3') == 1 ...
                     || strcmp(concStr, '4') == 1;
                % if it's a relevant concentration, find the pressure
                    pressStart = concEnd + 3;
                    pressureStr = file.name(pressStart:end-4);

                    if strcmp(pressureStr, '1') == 1 || ...
                       strcmp(pressureStr, '05') == 1 %|| ...
                     %  strcmp(pressureStr, '2') == 1         
                    % if it's a relevant pressure, load the file
                        % save the file name
                        dataA = fopen(strcat('D:\Cabs_Absorption Coefficients_h2o\', ...
                                      file.name));
                        dataB = textscan(dataA, '%f %f');

                        fileNameA = file.name(1:end-3); % cut off the txt
                        fileName = strcat(fileNameA, 'mat');

                        % store the cabs values
                        rawDataCabs = dataB{:, 2};

                        % generate the wavenumbers
                        rawDataWavenumbers = transpose(0:0.005:0.005*size(dataB{1,1},1)-0.005);%25000
                        
                        rawDataCabs = cat(2, rawDataWavenumbers, ...
                                        rawDataCabs);
                        
                        % save these values as a .mat file, using the old file 
                        % name as a base.
                        save(strcat('C:\Users\Scott\Documents\MATLAB\CabsDatabaseFiles\', ...
                                    fileName), 'rawDataCabs')
                    else
                    end
                else
                end
            else
            end       
        else
        end
    else
    end
end            
            