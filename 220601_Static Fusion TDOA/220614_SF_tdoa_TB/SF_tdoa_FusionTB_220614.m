clear; clc;close all
%% Simulation setup 
simTime = 30;                               % Total simulation time
simInterval = 5;                            % Update rate
[scenario, rxpairs, receiverIDs] = CreateMultiTargetTDOAScenario(simTime,simInterval); % Create scenario

%% Measurement statistics
measNoise = 1e0; % 100 ns per receiver
nFalse = 1; % no of false alarms per receiver pair
Pd = 0.99; % Detection probability per receiver pair




% Display object
display = TDOATrackingDisplay(XLimits=[-200e3 200e3],...
                                           YLimits=[-200e3 200e3],...
                                           LogAccuracy=false,...
                                           Title="TDOA Source Position Estimation");

% Create a GNN tracker
% tracker = trackerGNN(@initIMMFilter,...
%                      AssignmentThreshold=100);
tracker = trackerGNN( ...
    'FilterInitializationFcn',@initIMMFilter,...
    'MaxNumTracks', 20, ...
    'MaxNumSensors', 4, ...
    'AssignmentThreshold',100, ...
    'TrackLogic', 'Score', ...
    'DetectionProbability', 0.9, 'FalseAlarmRate', 1e-6, ...
    'Volume', 1e9, 'Beta', 1e-14);
% Define fuser. 
tdoaFuser = staticDetectionFuser(MaxNumSensors=receiverIDs(end)-1,...
    MeasurementFormat='custom',...
    MeasurementFcn=@MeasureTDOA,...
    MeasurementFusionFcn=@TDOA2Pos,...
    DetectionProbability=Pd,...
    FalseAlarmRate=1e-8,...
    UseParallel=false);

while advance(scenario)
    % Current elapsed time
    time = scenario.SimulationTime;

    % Simulate TDOA detections with false alarms and missed detections
    tdoaDets = SimulateTDOA(scenario, rxpairs, measNoise, false, Pd, nFalse);

    % Fuse TDOA detections to esimate position detections of unidentified
    % and unknown number of emitters
    posDets = tdoaFuser(tdoaDets);
    
    % Update tracker with position detections
    if ~isempty(posDets)
        tracks = tracker(posDets, time);
    end

    % Update display
    display(scenario, receiverIDs, tdoaDets, posDets,tracks);
end