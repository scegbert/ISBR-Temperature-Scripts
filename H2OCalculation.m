function [yh2o, yco2] = H2OCalculation(intMeasH2O, gasTempFinal, wallConditions, particleConditions,...
    pressure, yh2o, pathLength, H2OWavenumbers, fuelType)

yh2oPrevious = 0; %yh2o; %previous iteration of yh2o

if strcmp(fuelType, 'nat gas' ) == 1
    
    HtoC = 4; % natural gas
    H2OtoC = 0;
    
elseif strcmp(fuelType, 'fine wood' ) == 1
    
    HtoC = (5.36/1.008)/(49.87/12.011); %fine wood
    H2OtoC = (5.83/(2*1.008+15.999))/(49.87/12.011);
    
elseif strcmp(fuelType, 'med wood' ) == 1
    
    HtoC = (5.4/1.008)/(49.85/12.011); %medium wood
    H2OtoC = (5.28/(2*1.008+15.999))/(49.85/12.011);

else %unknown fuel type
    fuelType = throwerror; %this will throw an error if executed
    
end

yco2 = yh2o / (0.5 * HtoC + H2OtoC);

sensT = 300;
sensEm = 1;

keepGuessing = 1;
iter = 0;
limit = 10;
error(1) = 1;

while abs(error(1)) >= 0.001 && iter < limit

    % create new model from RTE equation
    [~, intsConv, nuConv] = GaussianConvolution(gasTempFinal, yh2o,yco2,...
            particleConditions(1), particleConditions(2), wallConditions(1),...
            wallConditions(2), sensT, sensEm, pressure, pathLength,...
            H2OWavenumbers(1), H2OWavenumbers(2));
    
    lambdaConv = 10^4 ./ nuConv;
    intConv(1) = trapz(flipud(lambdaConv),flipud(intsConv));

    %evaluate convergence
    error(2) = error(1); %slide previous iteration into 2nd array position
    error(1) = (intConv(1) - intMeasH2O) / intMeasH2O; %error between model and meas
    
    % guess another Yh2o
    if yh2oPrevious == 0
        % for the first iteration
        yh2oPrevious = yh2o;
        yh2o = yh2oPrevious * (1 - error(1));
    else
        % interpolate to 0
        yh2oInterpolated = (0-error(2)) / (error(1) - error(2)) * (yh2o - yh2oPrevious) + yh2oPrevious;
        yh2oPrevious = yh2o; %slide previous iteration into 2nd array position
        yh2o = yh2oInterpolated;
    end
     
    if yh2o >= .4 % error mitigation by bounding the potential results
        yh2o = .4;
        iter = iter + limit/2; %gives the algorithm one more chance to course correct
    elseif yh2o <= 0
        yh2o = 0;
        iter = iter + limit/2;
    end
            
    yco2 = yh2o / (0.5 * HtoC + H2OtoC);
    iter = iter + 1;
    intConv(2) = intConv(1);
  
end

end