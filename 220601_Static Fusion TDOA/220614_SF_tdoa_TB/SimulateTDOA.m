function detections = SimulateTDOA(scenario, rxpairs, measNoise, reportIdentity, Pd, nFalse)

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

nPairs = size(rxpairs,1);

% Find targets
platIDs = cellfun(@(x)x.PlatformID, scenario.Platforms);
targets = scenario.Platforms(~ismember(platIDs,rxpairs(:)));

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
sampleDetection = objectDetection(scenario.SimulationTime,0,'MeasurementParameters',measParams);

emissionSpeed = physconst('Lightspeed');
timeScale = 1e9;

detections = cell(0,1);

 for i = 1:size(rxpairs,1)
    sampleDetection.SensorIndex = i;
    r1Plat = rxpairs(i,1);
    r2Plat = rxpairs(i,2);
    plat1Pose = pose(scenario.Platforms{platIDs == r1Plat},'true');
    plat2Pose = pose(scenario.Platforms{platIDs == r2Plat},'true');
    interD = norm(plat2Pose.Position - plat1Pose.Position);
    maxTDOA = interD/emissionSpeed*timeScale;
    % True detections
    tdoa = zeros(numel(targets),1);
    identity = zeros(numel(targets),1);
    isDetected = rand(numel(targets),1) < Pd;
    for j = 1:numel(targets)
        targetPose = pose(targets{j},'true');
        r1 = norm(targetPose.Position - plat1Pose.Position);
        r2 = norm(targetPose.Position - plat2Pose.Position);
        trueTDOA = (r1 - r2)/emissionSpeed;
        reportedTDOA = trueTDOA*timeScale + sqrt(measNoise)*randn;
        if abs(reportedTDOA) > maxTDOA
            reportedTDOA = sign(reportedTDOA)*maxTDOA;
        end
        tdoa(j) = reportedTDOA;
        identity(j) = targets{j}.PlatformID;
    end
    tdoa = tdoa(isDetected);
    identity = identity(isDetected);

    % False detections
    tdoaFalse = maxTDOA*(-1 + 2*rand(nFalse(i),1));
    identityFalse = -1*ones(nFalse(i),1);
    identity = [identity;identityFalse]; %#ok<AGROW> 
    tdoa = [tdoa;tdoaFalse]; %#ok<AGROW> 

    % Fill detections
    thisTDOADetections = repmat({sampleDetection},numel(tdoa),1);
    for j = 1:numel(thisTDOADetections)
        thisTDOADetections{j}.Measurement = tdoa(j); 
        thisTDOADetections{j}.MeasurementNoise = measNoise;
        thisTDOADetections{j}.MeasurementParameters(1).OriginPosition = plat1Pose.Position(:);
        thisTDOADetections{j}.MeasurementParameters(2).OriginPosition = plat2Pose.Position(:);
        if reportIdentity
            thisTDOADetections{j}.ObjectClassID = identity(j);
        end
    end

    % From all receiver pairs
    detections = [detections;thisTDOADetections]; %#ok<AGROW> 
end

end