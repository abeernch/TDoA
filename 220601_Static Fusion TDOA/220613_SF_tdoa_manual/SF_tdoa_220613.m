clear;clc;close all
%%
sim_time = 32
time = 0;
numReceivers = 4;
theta = linspace(-pi,pi,numReceivers+1);
r = 5000;
xReceiver = r*cos(theta(1:end-1));
yReceiver = r*sin(theta(1:end-1));

% The algorithm is suitable for 3-D workflows. Each
% receiver must be at different height to observe/estimate the z of
% the object. 
zReceiver = zeros(1,numReceivers);
receiverpos = [xReceiver; yReceiver; zReceiver];

%% Measurement statistics
measNoise = 1e4; % 100 ns per receiver
nFalse = 2; % 2 false alarms per receiver pair
Pd = 0.95; % Detection probability per receiver pair
%% Add Target(s)
numTgts = 4;
xTarget = -5000 + 10000*rand(1,numTgts);
yTarget = -5000 + 10000*rand(1,numTgts);
vxTarget = 50*randn(1,numTgts);
vyTarget = 50*randn(1,numTgts);
zTarget = zeros(1,numTgts);
vzTarget = zeros(1,numTgts);
emitterpos = [xTarget; yTarget; zTarget];

%% Pairing sensors 
% PlatformIDs of TDOA calculation pairs. Each row represents the TDOA pair
% [1 3] means a TDOA is calculated between 1 and 3 with 3 as the reference
% receiver.
tdoaPairs = (1:(numReceivers-1))';
tdoaPairs(:,2) = numReceivers;

% IDs of all receivers
receiverIDs = 1:numReceivers;
reportIdentity = 0;

%% Define Static Fusion Algo object in MATLAB
tdoaFuser = staticDetectionFuser(MaxNumSensors=numReceivers-1,...
    MeasurementFormat='custom',...
    MeasurementFcn=@MeasureTDOA,...
    MeasurementFusionFcn=@TDOA2Pos,...
    DetectionProbability=Pd,...
    FalseAlarmRate=1e-8,...
    UseParallel=false);

%% Create Tracker 
tracker = trackerGNN(FilterInitializationFcn=@helperInitHighSpeedKF,...
                     AssignmentThreshold=100);
tracker.ConfirmationThreshold = [4 6];

%% initialize plot parameters
x = []; y = [];a= axes;xlim([-10e3 10e3]);ylim([-10e3 10e3])
scatter(receiverpos(1,:),receiverpos(2,:));
%% Initiate simloop
while time < sim_time + 0.1
    tdoaDets = SimulateTDOA(time,receiverpos,emitterpos,tdoaPairs, measNoise, false, Pd, nFalse,reportIdentity);
    
    % Fuse TDOA detections to esimate position detections of unidentified
    % and unknown number of emitters
    posDets = tdoaFuser(tdoaDets);
    
    % Update tracker with position detections
    if ~isempty(posDets)
        tracks = tracker(posDets, time);
    end
    time = time+1;
    plotTDOADetection(tdoaDets)
end 
