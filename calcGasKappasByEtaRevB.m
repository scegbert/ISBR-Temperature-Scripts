function [kappa] = calcGasKappasByEtaRevB(gasName, pressure, gasT, conc, ...
                    wavenumberStart, wavenumberStop, wavenumbersFTIR, ...
                    rowLowAFTIR, rowHighEFTIR)
                    %testNum,
% Author: John Tobiasson
% Date: 04/07/2016
% Purpose: This function will caclculate an array of gas emissivities for a
% given set of conditions, from John Person's Cabs data, created with Dr. 
% Solovjov. The resulting array will be passed to another program. The
% Cabs used will be created from interpolating between temperature and
% concentration files.

% Necesary Inputs:
% pressure = total pressure in atm
% gasT = gas temperature in K
% conc = water concentration, mole fraction
% wavenumberStart = beginning wavenumber in region of interest, 1/cm
% wavenumberStop = ending wavenumber in region of interest, 1/cm

% Optional Inputs:
% wavenumbersFTIR = set of wavenumbers to pull information for, if
    % spacing other than 0.005 1/cm is used.
% rowLowFTIR = row in wavenumbersFTIR to start with
% rowHighFTIR = row in wavenumbersFTIR to end with

% RevB note: This revision will include interpolation on pressure, and will
% include an allowance for not inputing FTIR wavenumbers. If no FTIR
% wavenumbers are included, the code will assume that standard wavenumbers
% will be used (5185:0.005:5715). The spacing of 0.005  is from Pearson's 
% database. 5185 is the cm^-1 where my gas temperature processing starts,
% and 5715 is where it ends, spectrally.

% RevC note: This revision changes from bringing in Cabs from .txt files,
% and loads them from .mat files. This changes from taking 6-25s, depending
% on the load method used, to only taking on the order of 1s , sometimes
% less.

% define constants
%pressure = 0.84445041; % atm, based on BYU elevation of  aprox 4650 ft, 12.41 psia
Ru = 0.8205; % universal gas constant, cm^2 * m / mol / K
Na = 6.02214085774 * (10^23); % avogadro's number, molecules / mol

calcConc = conc; % concentration to be used.
calcPressure = pressure; % pressure to be used.

if strcmp(gasName, 'h2o') == 1
    % location of the H2O absorption cross section database (Pearson)
    fileLoc = 'F:\relevantH2OCabsDatabaseFiles\';

else
    % location of the CO2 absorption cross section database
    fileLoc = 'F:\relevantCO2CabsDatabaseFiles\';
    
end

databaseLoc = 'E:\Cabs_Absorption Coefficients_'; % This is where the .txt
    % database files are located, the ones that Pearson created. My files
    % are actually located in a place like this: 'E:\Cabs_Absorption
    % Coefficients_H2O\' The H2O absorption cross section files are then
    % stored in that directory. I only take things this far here because I
    % have separate directories for H2O and CO2, and both will be needed. I
    % handle that later in the code.
    % This is necessary for attempting to load in, and then convert to
    % binary, files that are still only found in a .txt version that
    % Pearson generated.

checkForEtaInputs = exist('wavenumberStart','var'); % check for the existence
                    % of a wavenumberStart variable
if checkForEtaInputs == 1
    checkForEtaInputs = exist('wavenumberStop','var'); % check for the existence
                    % of a wavenumberStop variable
    if checkForEtaInputs == 1 % If these inputs were given, use them.
        etaLow = wavenumberStart;
        etaHigh = wavenumberStop;
    else
    end
else
    etaLow = 5185; % Beginning of H2O property processing region of interest
    etaHigh = 5715; % Ending of H2O property processing region of interest
end

% determine files to use for linear interpolation. This allows the function
% to interpolate based on pressure, tempearture, and concentration (molar).
% These values are based on what is available in Pearson's database.
    % determine high and low temperatures
    tLow = 100 .* floor(gasT ./ 100);
    tHigh = 100 .* ceil(gasT ./ 100);
     
    % determine high and low concentrations
    if strcmp(gasName, 'h2o')
        concDet = 10 .* conc;

        if concDet == 0

            concNums = 0;
            concStrs = {'0'};

        elseif concDet == 0.5

            concNums = 0.05;
            concStrs = {'05'};

        elseif concDet == 1

            concNums = 0.1;
            concStrs = {'1'};

        elseif concDet == 2

            concNums = 0.2;
            concStrs = {'2'};

        elseif concDet == 3

            concNums = 0.3;
            concStrs = {'3'};

        elseif concDet == 4

            concNums = 0.4;
            concStrs = {'4'};      

        elseif concDet == 6

            concNums = 0.6;
            concStrs = {'6'};

        elseif concDet == 8

            concNums = 0.8;
            concStrs = {'8'};

        elseif concDet == 10

            concNums = 1.0;
            concStrs = {'10'};

        else

            if concDet < 0.5

                concNums = [0; 0.05];
                concStrs = {'0'; '05'};

            elseif concDet < 1

                concNums = [0.05; 0.1];
                concStrs = {'05'; '1'};

            elseif concDet < 2

                concNums = [0.1; 0.2];
                concStrs = {'1'; '2'};

            elseif concDet < 3

                concNums = [0.2; 0.3];
                concStrs = {'2'; '3'};

            elseif concDet < 4

                concNums = [0.3; 0.4];
                concStrs = {'3'; '4'};

            elseif concDet < 6

                concNums = [0.4; 0.6];
                concStrs = {'4'; '6'};

            elseif concDet < 8

                concNums = [0.6; 0.8];
                concStrs = {'6'; '8'};

            elseif concDet < 10

                concNums = [0.8; 1.0];
                concStrs = {'8'; '10'};

            end

        end
    
    else        
        % gas is CO2, and all concentrations are 0 for looking up the
        % files. This is important because Pearson's database does not
        % provide different values for the different concentrations of CO2,
        % as he found it to be a negligible parameter in determining CO2
        % absorption cross sections.
        concNums = 0;
        concStrs = {'0'};
    end
    if strcmp(gasName, 'h2o') == 1
        
        if pressure == 0.1

            pressNums = 0.1;
            pressStrs = {'01'};

        elseif pressure == 0.25

            pressNums = 0.25;
            pressStrs = {'025'};

        elseif pressure == 0.5

            pressNums = 0.5;
            pressStrs = {'05'};

        elseif pressure == 1

            pressNums = 1;
            pressStrs = {'1'};

        elseif pressure == 2

            pressNums = 2;
            pressStrs = {'2'};

        elseif pressure == 4

            pressNums = 4;
            pressStrs = {'4'};      

        elseif pressure == 8

            pressNums = 8;
            pressStrs = {'8'};

        elseif pressure == 15

            pressNums = 15;
            pressStrs = {'15'};

        elseif pressure == 30

            pressNums = 30;
            pressStrs = {'30'};

        elseif pressure == 50

            pressNums = 50;
            pressStrs = {'50'};

        else

            if pressure < 0.1

                pressNums = [0; 0.1];
                pressStrs = {'0'; '01'};

            elseif pressure < 0.25

                pressNums = [0; 0.25];
                pressStrs = {'0'; '025'};

            elseif pressure < 0.5

                pressNums = [0; 0.5];
                pressStrs = {'0'; '05'};

            elseif pressure < 1

                pressNums = [0.5; 1];
                pressStrs = {'05'; '1'};

            elseif pressure < 2

                pressNums = [1; 2];
                pressStrs = {'1'; '2'};

            elseif pressure < 4

                pressNums = [2; 4];
                pressStrs = {'2'; '4'};

            elseif pressure < 8

                pressNums = [4; 8];
                pressStrs = {'4'; '8'};

            elseif pressure < 15

                pressNums = [8; 15];
                pressStrs = {'8'; '15'};


            elseif pressure < 30

                pressNums = [15; 30];
                pressStrs = {'15'; '30'};

            else

                pressNums = [30; 50];
                pressStrs = {'30'; '50'};

            end

        end       
    else
        pressNums = 1;
        pressStrs = '1';
    end

% From here on, this is one massive linear interpolation scheme.    
    
% load Cabs files
if size(pressNums, 1) == 1

    if size(concNums, 1) == 1
    
        if tHigh == tLow
        
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc, gasName,'\', ...
                    filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');

            end
             
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = size(rawDataCabs1, 1);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
            
            % dimension and fill array
            clear('rawDataCabs');
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
            rawDataCabs(:, 2) = rawDataCabs1(rowCabsLow:rowCabsHigh, 2);
        
        else
        
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);

            end
    
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs2 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs2 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs2(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs2(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs2');
                %rawDataCabs2 = rawDataCabs2(:,2);

            end
    
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs2, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
    
            % dimension array
            clear('rawDataCabs');
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
    
            % interpolate between Cabs by temperature
            r = 1;
    
            for i = rowCabsLow:1:rowCabsHigh;
        
                rawDataCabs(r, 2) = (gasT - tLow) ./ (tHigh - tLow) .* ...
                    (rawDataCabs2(i, 2) - rawDataCabs1(i, 2)) + rawDataCabs1(i, 2);
        
                r = r + 1;
    
            end
        
        end
    
    else
    
    % interpolate by concentration, then by temperature
        if tHigh == tLow
        
            % only need to interpolate by concentration
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});           
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(2), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs3 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs3 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs3(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs3(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs3');
                %rawDataCabs3 = rawDataCabs3(:,2);

            end
        
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs3, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
    
            % dimension arrays
            clear('rawDataCabs');
            rawDataCabs = zeros(rowCabsHigh - rowCabsLow + 1, 2);
    
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
    
            q = 1;
    
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabs(q, 2) = (conc - concNums(1)) ./ (concNums(2) - ...
                    concNums(1)) .* (rawDataCabs3(j, 2) - rawDataCabs1(j, 2)) + ...
                    rawDataCabs1(j, 2);
        
                q = q + 1;
        
            end
        
        else
            % interpolate by concentration, then temperature
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs2 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs2 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs2(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs2(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs2');
                %rawDataCabs2 = rawDataCabs2(:,2);
                
            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(2), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs3 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs3 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs3(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs3(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs3');
                %rawDataCabs3 = rawDataCabs3(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(2), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs4 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs4 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs4(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs4(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs4');
                %rawDataCabs4 = rawDataCabs4(:,2);

            end
    
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs2, 1), ...
                                size(rawDataCabs3, 1), size(rawDataCabs4, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
    
            % dimension arrays
            clear('rawDataCabs');
            rawDataCabs = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt1 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt2 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
    
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
    
            q = 1;
    
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabsInt1(q, 2) = (conc - concNums(1)) ./ (concNums(2) - ...
                    concNums(1)) .* (rawDataCabs3(j, 2) - rawDataCabs1(j, 2)) + ...
                    rawDataCabs1(j, 2);
        
                rawDataCabsInt2(q, 2) = (conc - concNums(1)) ./ (concNums(2) - ...
                    concNums(1)) .* ...
                    (rawDataCabs4(j, 2) - rawDataCabs2(j, 2)) + rawDataCabs2(j, 2);
        
                q = q + 1;
        
            end
    
            % interpolate by temperature
            for k = 1:1:size(rawDataCabsInt1, 1);
        
                rawDataCabs(k, 2) = (gasT - tLow) ./ (tHigh - tLow) .* ...
                    (rawDataCabsInt2(k, 2) - rawDataCabsInt1(k, 2)) + ...
                    rawDataCabsInt1(k, 2);
    
            end
        end
    end
    
else
    
    if size(concNums, 1) == 1
    
        if tHigh == tLow
            % interpolate by pressure only
            filenameCabs = strcat({fileLoc}, {'Cabs'},{' '},{gasName},{' '}, ...
                {num2str(tLow)}, {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});           
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);
                
            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1}, 'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs5 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs5 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs5(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs5(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs5');
                %rawDataCabs5 = rawDataCabs5(:,2);

            end
            
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs5, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
    
            % dimension array
            clear('rawDataCabs');
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
            
            % interpolate by pressure
            q = 1;
            for k = rowCabsLow:1:rowCabsHigh
                rawDataCabs(q, 2) = (pressure - pressNums(1)) ./ ...
                                (pressNums(2) - pressNums(1)) .* ...
                    (rawDataCabs5(k, 2) - rawDataCabs1(k, 2)) + ...
                    rawDataCabs1(k, 2);
                
                q = q + 1;
                
            end
    
            % fill array
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
            rawDataCabs(:, 2) = rawDataCabs1(rowCabsLow:rowCabsHigh, 2);
        
        else
            % interpolate by pressure, then by temperature
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs2 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs2 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs2(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs2(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs2');
                %rawDataCabs2 = rawDataCabs2(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs5 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs5 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs5(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs5(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs5');
                %rawDataCabs5 = rawDataCabs5(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(1), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs6 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs6 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs6(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs6(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs6');
                %rawDataCabs6 = rawDataCabs6(:,2);
                
            end
    
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs2, 1), ...
                                size(rawDataCabs5, 1), size(rawDataCabs6, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
            % dimension arrays
            clear('rawDataCabs');
            rawDataCabs = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt1 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt2 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
    
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);

            % interpolate by pressure
            q = 1;
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabsInt1(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* (rawDataCabs5(j, 2) - rawDataCabs1(j, 2)) + ...
                    rawDataCabs1(j, 2);
        
                rawDataCabsInt2(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* ...
                    (rawDataCabs6(j, 2) - rawDataCabs2(j, 2)) + rawDataCabs2(j, 2);
        
                q = q + 1;
        
            end
   
            % interpolate between Cabs by temperature
            for k = 1:1:size(rawDataCabsInt1, 1);
        
                rawDataCabs(k, 2) = (gasT - tLow) ./ (tHigh - tLow) .* ...
                    (rawDataCabsInt2(k, 2) - rawDataCabsInt1(k, 2)) + ...
                    rawDataCabsInt1(k, 2);
    
            end
        
        end
    
    else
    
        if tHigh == tLow
        
            % interpolate by pressure, then by concentration
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(2), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs3 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs3 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs3(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs3(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs3');
                %rawDataCabs3 = rawDataCabs3(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs5 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs5 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs5(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs5(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs5');
                %rawDataCabs5 = rawDataCabs5(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(2), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs7 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs7 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs7(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs7(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs7');
                %rawDataCabs7 = rawDataCabs7(:,2);

            end
            
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs3, 1), ...
                                size(rawDataCabs5, 1), size(rawDataCabs7, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
    
            % dimension arrays
            clear('rawDataCabs');
            rawDataCabs = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt1 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt2 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
    
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
    
            % interpolate by pressure
            q = 1;
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabsInt1(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* (rawDataCabs5(j, 2) - rawDataCabs1(j, 2)) + ...
                    rawDataCabs1(j, 2);
        
                rawDataCabsInt2(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* ...
                    (rawDataCabs7(j, 2) - rawDataCabs3(j, 2)) + rawDataCabs3(j, 2);
        
                q = q + 1;
        
            end
            
            % interpolate by concentration
            q = 1;
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabs(q, 2) = (conc - concNums(1)) ./ (concNums(2) - ...
                    concNums(1)) .* (rawDataCabsInt2(q, 2) - rawDataCabsInt1(q, 2)) + ...
                    rawDataCabsInt1(q, 2);
        
                q = q + 1;
        
            end
        
        else
            % interpolate by pressure, then by concnetration, then by
            % temperature
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs1 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs1 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs1(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs1(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs1');
                %rawDataCabs1 = rawDataCabs1(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(1), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs2 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs2 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs2(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs2(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs2');
                %rawDataCabs2 = rawDataCabs2(:,2);

            end
    
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(2), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs3 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs3 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs3(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs3(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs3');
                %rawDataCabs3 = rawDataCabs3(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(2), {' P'}, pressStrs(1),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs4 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs4 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs4(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs4(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs4');
                %rawDataCabs4 = rawDataCabs4(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(1), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs5 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs5 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs5(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs5(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs5');
                %rawDataCabs5 = rawDataCabs5(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(1), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs6 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs6 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs6(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs6(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs6');
                %rawDataCabs6 = rawDataCabs6(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tLow), {'K Y'}, concStrs(2), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs7 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs7 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs7(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs7(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs7');
                %rawDataCabs7 = rawDataCabs7(:,2);

            end
            
            filenameCabs = strcat({fileLoc}, {'Cabs '},{gasName},{' '}, ...
                num2str(tHigh), {'K Y'}, concStrs(2), {' P'}, pressStrs(2),...
                {'.mat'});
            % check if the file exists. If not, try to load the .txt
            % version- if successful here, save it as a .mat
            if exist(filenameCabs{1},'file') == 2;    
                
                rawDataCabs = load(filenameCabs{1});
                fieldName = fields(rawDataCabs);
                rawDataCabs8 = rawDataCabs.(fieldName{1});
                
            else
                
                loadFilename = strcat(databaseLoc,...
                    gasName,'\',filenameCabs{1}((size(fileLoc,2)+1):(end-3)),'txt');
                fileID = fopen(loadFilename);
                rawDataCabs = textscan(fileID,'%s %s');
                rawDataCabs8 = zeros(size(rawDataCabs{1,1},1), 2);
                rawDataCabs8(:,1) = str2double(rawDataCabs{1, 1});
                rawDataCabs8(:,2) = str2double(rawDataCabs{1, 2});
                save(filenameCabs{1},'rawDataCabs8');
                %rawDataCabs8 = rawDataCabs8(:,2);

            end
            
            % find necessary rows for Cabs files
            [~, rowCabsLow] = min(abs(rawDataCabs1(:, 1) - etaLow));            
            if etaHigh == 0
                rowCabsHigh = min([size(rawDataCabs1, 1), size(rawDataCabs2, 1), ...
                                size(rawDataCabs3, 1), size(rawDataCabs4, 1),...
                                size(rawDataCabs5, 1), size(rawDataCabs6, 1), ...
                                size(rawDataCabs7, 1), size(rawDataCabs8, 1)]);
            else
                [~, rowCabsHigh] = min(abs(rawDataCabs1(:, 1) - etaHigh));
            end
    
            % dimension arrays
            clear('rawDataCabs');
            rawDataCabs = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt1 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsInt2 = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsIntA = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsIntB = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsIntC = zeros(rowCabsHigh - rowCabsLow + 1, 2);
            rawDataCabsIntD = zeros(rowCabsHigh - rowCabsLow + 1, 2);
    
            rawDataCabs(:, 1) = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
    
            % interpolate by pressure
            q = 1;
    
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabsIntA(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* (rawDataCabs5(j, 2) - rawDataCabs1(j, 2)) + ...
                    rawDataCabs1(j, 2);
        
                rawDataCabsIntB(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* ...
                    (rawDataCabs6(j, 2) - rawDataCabs2(j, 2)) + rawDataCabs2(j, 2);
        
                rawDataCabsIntC(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* (rawDataCabs7(j, 2) - rawDataCabs3(j, 2)) + ...
                    rawDataCabs3(j, 2);
        
                rawDataCabsIntD(q, 2) = (pressure - pressNums(1)) ./ (pressNums(2) - ...
                    pressNums(1)) .* ...
                    (rawDataCabs8(j, 2) - rawDataCabs4(j, 2)) + rawDataCabs4(j, 2);
                
                q = q + 1;
        
            end
            
            % interpolate by concentration
            q = 1;
    
            for j= rowCabsLow:1:rowCabsHigh;
        
                rawDataCabsInt1(q, 2) = (conc - concNums(1)) ./ (concNums(2) - ...
                    concNums(1)) .* (rawDataCabsIntC(q, 2) - rawDataCabsIntA(q, 2)) + ...
                    rawDataCabsIntA(q, 2);
        
                rawDataCabsInt2(q, 2) = (conc - concNums(1)) ./ (concNums(2) - ...
                    concNums(1)) .* ...
                    (rawDataCabsIntD(q, 2) - rawDataCabsIntB(q, 2)) + rawDataCabsIntB(q, 2);
        
                q = q + 1;
        
            end                
            
            % interpolate by temperature
            for k = 1:1:size(rawDataCabsInt1, 1);
        
                rawDataCabs(k, 2) = (gasT - tLow) ./ (tHigh - tLow) .* ...
                    (rawDataCabsInt2(k, 2) - rawDataCabsInt1(k, 2)) + ...
                    rawDataCabsInt1(k, 2);
    
            end
        end
    end
    
end

% This ends the massive interpolation scheme.

% save Cabs
Cabs = rawDataCabs(:, 2);

% create wavenumbers arrays
wavenumbersCabs = rawDataCabs1(rowCabsLow:rowCabsHigh, 1);
    
% calculate gas emissivities for each wavenumber

% check if FTIR wavenumbers were passed to the function. This indicates
% whether or not this is simply for generating a model, or to be used in 
% the processing an actual measurement, which will have different spacing.
% For simply generating a model, this is not necessary, not is it used.
checkForArray = exist('wavenumbersFTIR', 'var'); % wavenumbersFTIR would be
    % an array of wavenumbers coming from a measurement with the FTIR.

if checkForArray == 1
    % dimension emissivities array
    gasEms = zeros(rowHighEFTIR - rowLowAFTIR + 1, 1);

    v = 1;

    for m = rowLowAFTIR:1:rowHighEFTIR
    
        % find wavenumber row in Cabs data
        [~, rowCabs] = min(abs(wavenumbersCabs - wavenumbersFTIR(m)));
    
        % calculate kappa
        kappa(:,2) = Cabs(rowCabs) .* ((calcConc * calcPressure) ./ (Ru .* gasT)) * Na;
    
    end
    
    wavenumbersToWrite = wavenumbersFTIR;
    
    kappa(:, 1) = wavenumbersToWrite;
    
else    
    
    kappa(:, 2) = Cabs .* ((calcConc * calcPressure) ./ (Ru .* gasT)) * Na;
    
    wavenumbersToWrite = wavenumbersCabs;
    
    kappa(:, 1) = wavenumbersToWrite;
    
end

%xlswrite('testGasEmsWithPsTsYs.xlsx', cat(2, wavenumbersToWrite, gasEms), testNum);

end