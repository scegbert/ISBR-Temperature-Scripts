function [intsModel, intsConv, nuFinal] = GaussianConvolution(...
                                    gasT, gasConcH2O, gasConcCO2, partT, partKappa,...
                                    wallT, wallEm, sensT, sensEm, pressure,...
                                    pathLength, startWavenumber, stopWavenumber)

    %creates spectral model from RTE equation (function has same inputs as RTE)
    %then takes the convolution of that spectral modelm, returning both
                                
    [intsModel, nuModel] = OneD_RTE_GAP(gasT, gasConcH2O, gasConcCO2, partT, partKappa,...
                                    wallT, wallEm, sensT, sensEm, pressure,...
                                    pathLength, startWavenumber-1, stopWavenumber+1);
                                    %expands wavenumber range prior to convolution
    
    numberOfWeights = 70; % must be an even number for the weighting to work properly
    sigma = 1;
    numSigmasFromCenter = 1;
    stepsize = 2 * numSigmasFromCenter * sigma / numberOfWeights;

    x = transpose(-numSigmasFromCenter*sigma:stepsize:numSigmasFromCenter*sigma);   
    norms = (1 / ((2 * (sigma^2) * pi)^0.5)) .* exp(-(x.^2)./(2*(sigma^2)));
    w = norms ./ sum(norms);    

    % loop through model spectrum, weighting sections and averaging them 
        % find where to start the weighting and where to end it
        weightSize = size(w, 1);
        startWeighting = 1 + floor(weightSize / 2);
        stopWeighting = size(intsModel, 1) - floor(weightSize / 2);

        intsConv = zeros((size(intsModel, 1) - 2 * floor(weightSize / 2)), 1);
        nuConv = zeros(size(intsConv, 1), 1);

        newIter = 1;

    for j = startWeighting:1:stopWeighting

        nuConv(newIter, 1) = nuModel(j);
        intsConv(newIter, 1) = sum(w .* intsModel(j-floor(weightSize / 2): ...
                        j+floor(weightSize / 2)), 1);

        newIter = newIter + 1;

    end
    
    %chop model ints to match measured ints
    [~, removeLow] = min(abs(nuModel-startWavenumber));
    [~, removeHigh] = min(abs(nuModel-stopWavenumber));    
    intsModel = intsModel(removeLow:removeHigh); 

    %chop model wavenumbers to match the model ints
    nuFinal = nuModel(removeLow:removeHigh);
    
    %chop convoluted ints to match modeled ints
    [~, removeLow] = min(abs(nuConv-startWavenumber));
    [~, removeHigh] = min(abs(nuConv-stopWavenumber));    
    intsConv = intsConv(removeLow:removeHigh); 
    
end