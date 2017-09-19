function [intsModel, nuModel] = OneD_RTE_MOD( ...
    tempGas, yh2o, yco2, tempPart, kappaPart,...
    tempWall, emrho, pressure, pathLength,...
    startWavenumber, stopWavenumber)

%% set up spectral properties and heat transfer boundary conditions

% wavenumbers of interest
nuModel = transpose(startWavenumber:0.005:stopWavenumber);

tauSens = .88; %transmisivity of sapphire in IR spectrum, Thorlabs
emSens = (1-tauSens) / 2; %absorbed and emitted portion, assumed = rho
rhoSens = (1-tauSens) / 2; %reflected portion, assumed = em
tempSens = 300; %Temperature of the sapphire window

tempTarget = 300; %temperature of the cold target

emWall = .7; %emissivity of the BFR wall;

%emrho = emWall * rhoTarget = the value measured uding 2 color pyrometry
rhoTarget = emrho / emWall; %attentuation of wall emission due to reflection
emTarget = 1 - rhoTarget; %opaque wall

% black body intensities
intsSens = PlanckIntensities(tempSens, nuModel) * emSens;
intsWall = PlanckIntensities(tempWall, nuModel) * emWall;
intsTarget = PlanckIntensities(tempTarget, nuModel) * emTarget;

%% verify and calculate gas properties
%verify that input values are either constant or the same length
inputLengths = unique([length(tempGas), length(yh2o), length(yco2), length(tempPart), length(kappaPart)]);

if length(inputLengths) > 2 %should be length 1 (constant) with max of 1 other length (varying w/ PL)
    inputLengths = throwerror; %this will throw an error if multiple arrays are different lengths
end

ni = max(inputLengths);
dx = pathLength/ni;

% make all constant input properties the same length (ni)
if length(tempGas) == 1
    tempGas = ones(1, ni) * tempGas;
end

if length(yh2o) == 1
    yh2o = ones(1, ni) * yh2o;
end

if length(yco2) == 1
    yco2 = ones(1, ni) * yco2;
end

if length(tempPart) == 1
    tempPart = ones(1, ni) * tempPart;
end

if length(kappaPart) == 1
    kappaPart = ones(1, ni) * kappaPart; %assuming gray particles
end

for j = 1:1:ni % gas properties for each cell
    %if tempGas yh2o are the same as previous, use previous values
    if j > 1 && tempGas(j) == tempGas(j-1) && yh2o(j) == yh2o(j-1)
        kappaH2O(:, j) = kappaH2O(:, j - 1);
        
    else
        [kappaH2O(:, j), ~ ] = calcGasKappasByEtaRevB('h2o', pressure, tempGas(j),...
            yh2o(j), startWavenumber, stopWavenumber);
        
    end
    
    %if tempGas yco2 are the same as previous, use previous values
    if j > 1 && tempGas(j) == tempGas(j-1) && yco2(j) == yco2(j-1)
        kappaCO2(:, j) = kappaCO2(:, j - 1);
        
    else
        
        [kappaCO2(:, j), ~ ] = calcGasKappasByEtaRevB('co2', pressure,...
            tempGas(j), yco2(j), startWavenumber, stopWavenumber);
        
        %% abskgCO2iter(end+1,size(abskgH2Oiter,2),:) = 0; %add zeros to make them same length - send to kappa function
        
    end
    
    %combinations of kappas
    kappaGas(:,j) = kappaH2O(:, j) + kappaCO2(:, j);
    kappaAll(:,j) = kappaH2O(:, j) + kappaCO2(:, j) + kappaPart(j);
    
    % emissive optical thicknesses
    tauGas(:,j) = kappaGas(j) * dx(j);
    tauPart(:,j) = kappaPart(j) * dx(j);
    
    % absorptive optical thickness
    tauAbs(:,j) = kappaAll(j) * dx(j);
    
    % black body intensities
    intsPart(:,j) = PlanckIntensities(tempPart(j), nuModel) .* kappaPart(j);
    intsGas(:,j) = PlanckIntensities(tempGas(j), nuModel) .* kappaGas(j);
    
end

%% Start sweep from sensor side (L, i = 0) toward the target (R, i = ni)

qinL = 0; % incident flux on the left (sensor, i = 0) face
dqinL = 1; % change in incident flux on the left (sensor, i = 0) face
dqinR = 1; % change in incident flux on the right (target, i = ni) face

while dqinL > 0.0001 || dqinR > 0.0001 %convergence criteria, both must be met
    
    for i = 1:1:ni %starting at the sensor, working out to wall (and target)
        
        if i == 1
            % entering first cell (i = 1), intensity = emitted + reflected off sensor
            IoutR(:,i) = intsSens + intsWall + rhoSens .* qinL / pi;
        else
            % intensity entering cell = intensity leaving previous
            IinR(:,i) = IoutR(:,i-1);
        end
        % intensity (moving in +i) generated/absorbed in each cell
        % (see Tobiasson thesis equations 5-1, -8, -9)
        IoutR(:,i) = IinR(:, i) .* exp(-tauAbs) + ... %absorption
            (1-exp(-tauAbs) ./  kappaAll) .* (intsGas + intsPart); %emission
        
    end
    qinR = IoutR(:, 1) * pi; % incident flux on the right (target, i = ni) face
    
    for i = ni:1:1 % starting at the target, working back to sensor
        
        if i == ni
            % entering first cell (i = ni),
            % intensity = emitted + reflected off wall and target
            IoutL(:,i) = intsTarget + rhoTarget .* qinR / pi;
        else
            % intensity entering cell = intensity leaving previous
            IinL(:,i) = IoutL(:,i+1);
        end
        % intensity (moving in -i) generated/absorbed in each cell
        % (see Tobiasson thesis equations 5-1, -8, -9)
        IoutL(:,i) = IinL(:, i) .* exp(-tauAbs) + ... %absorption
            (1-exp(-tauAbs) ./  kappaAll) .* (intsGas + intsPart); %emission
        
    end
    
    qinL = IoutL(:, 1) * pi; %incident flux at sensor
    
    % Check for convergence.
    dqinL = max((qinL - qinLprev) ./ qinLprev); % max change in incident flux
    dqinR = max((qinR - qinRprev) ./ qinRprev);
    
    qinRprev = qinR; %previuos value of incident flux at target
    qinLprev = qinL; %previous value of incident flux at sensor
    
end

IntsModel = (qinL / pi) * tauSens;

end
