function detections = SimulateTDOA(time,receiverpos,emitterpos,rxPairs, measNoise, false, Pd, nFalse,reportIdentity)
% This function simulates TDOA for each TDOA pair defined in the input,
% tdoaPairs. 
%
% Each property of the objectDetection object is defined as
% 
% SensorIndex - Unique identifier of each pair. 
% Time - Current elapsed simulation time (s)
% Measurement - TDOA measurement of target (ns)
% MeasurementNoise - uncertainty variance in TDOA measurement (ns^2)
% MeasurementParameters - A 2-by-1 struct defining the position of
% receivers in the scenario. Each struct has field "OriginPosition" which
% represents the position of the receivers. 

nPairs = size(rxPairs,1);
nEmitter = size(emitterpos,2);

% Find targets
% platIDs = cellfun(@(x)x.PlatformID, scenario.Platforms);
% targets = scenario.Platforms(~ismember(platIDs,rxPairs(:)));

if nargin < 3
    measNoise = 0;
end

if nargin < 4
    reportIdentity = false;
end

if nargin < 5
    Pd = 1;
end

if nargin < 6
    nFalse = 0;
end

% False alarms per pair
if isscalar(nFalse)
    nFalse = nFalse*ones(nPairs,1);
end

% Define sampleDetection
measParams = repmat(struct('OriginPosition',zeros(3,1)),2,1);
sampleDetection = objectDetection(time,0,'MeasurementParameters',measParams);

emissionSpeed = physconst('LightSpeed');
timeScale = 1e9;

detections = cell(0,1);

for i = 1:size(rxPairs,1)
    sampleDetection.SensorIndex = i;
    r1Plat = rxPairs(i,1);
    r2Plat = rxPairs(i,2);
    plat1Pos = receiverpos(:,r1Plat);
    plat2Pos = receiverpos(:,r2Plat);
    interD = norm(plat2Pos - plat1Pos);
    maxTDOA = interD/emissionSpeed*timeScale;
    % True detections
    tdoa = zeros(numel(nEmitter),1);
    identity = zeros(numel(nEmitter),1);
    isDetected = rand(numel(nEmitter),1) < Pd;
    for j = 1:numel(nEmitter)
        targetPos = emitterpos(:,j);
        r1 = norm(targetPos - plat1Pos);
        r2 = norm(targetPos - plat2Pos);
        trueTDOA = (r1 - r2)/emissionSpeed;
        reportedTDOA = trueTDOA*timeScale + sqrt(measNoise)*randn;
        if abs(reportedTDOA) > maxTDOA
            reportedTDOA = sign(reportedTDOA)*maxTDOA;
        end
        tdoa(j) = reportedTDOA;
        identity(j) = j;
    end
    tdoa = tdoa(isDetected);
    identity = identity(isDetected);

    % False detections
    tdoaFalse = maxTDOA*(-1 + 2*rand(nFalse(i),1));
    identityFalse = 2*ones(nFalse(i),1);
    identity = [identity;identityFalse]; %#ok<AGROW> 
    tdoa = [tdoa;tdoaFalse]; %#ok<AGROW> 

    % Fill detections
    thisTDOADetections = repmat({sampleDetection},numel(tdoa),1);
    for j = 1:numel(thisTDOADetections)
        thisTDOADetections{j}.Measurement = tdoa(j); 
        thisTDOADetections{j}.MeasurementNoise = measNoise;
        thisTDOADetections{j}.MeasurementParameters(1).OriginPosition = plat1Pos;
        thisTDOADetections{j}.MeasurementParameters(2).OriginPosition = plat2Pos;
        if reportIdentity
            thisTDOADetections{j}.ObjectClassID = identity(j);
        end
    end

    % From all receiver pairs
    detections = [detections;thisTDOADetections]; %#ok<AGROW> 
end

end