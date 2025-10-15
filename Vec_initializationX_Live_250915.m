%cd ('\Data\4 nodes')
clc; clear; close all;

addpath('D:\SC2\Projects\TDOA\MS17 TDoA Integrated Scripts\Matlab_testing\Data\4 nodes\')
addpath('D:\SC2\Projects\TDOA\MS17 TDoA Integrated Scripts\Matlab_testing\Functions')

filename = 'tdoa_mar15_4_7609EC_both_filters.mat';
load(filename)
% profile viewer
% profile on

tic;

%% Convert Sensor LAT,LON to local
% Node Positions in geographic coordinates
% Node 1: NASTP
% Node 2: Swana
% Node 3: I-14
% Node 4: Chak Shahzad
[lat_node{1},lon_node{1},alt_node{1}] = deal(33.609144,73.104425,510);
[lat_node{2},lon_node{2},alt_node{2}] = deal(33.5438077,73.1379735,439);
[lat_node{3},lon_node{3},alt_node{3}] = deal(33.6088889,72.9805555,567);
[lat_node{4},lon_node{4},alt_node{4}] = deal(33.674690,73.184789,536);

no_of_nodes = 4;
nSensors = size(lat_node,2);
rx = cell(1,nSensors);
SensorName = {'NASTP', 'Swana', 'I-14','ChakShahzad'};
sensorMask = true(1,no_of_nodes);
sensorInd = find(sensorMask);
% sensorMask = logical(round(rand(1,4))) ;
% sensorInd = find(sensorMask==1);

%% Initialize Sensor Position
for i = 1:nSensors
    [rx{i}.pos(1,1),rx{i}.pos(2,1),rx{i}.pos(3,1)] = ...
        latlon2local(lat_node{i},lon_node{i},alt_node{i},[lat_node{1},lon_node{1},alt_node{1}]);
end

%% Active sensors
% Select and initialize the active/enabled sensors based on received mask (4,3,2)

if sum(sensorMask) == 4
    sensor_pos = [rx{1}.pos rx{2}.pos rx{3}.pos rx{4}.pos];
    
    tdoaPairs(:,2) = (2:nSensors)'; tdoaPairs(:,1) = 1;

    syms x1 x2 x3 x4 xx y1 y2 y3 y4 yy z1 z2 z3 z4 zz c tdoa12 tdoa13 tdoa14
    eqn1 = sqrt((x1-xx)^2+(y1-yy)^2+(z1-zz)^2) - sqrt((x2-xx)^2+(y2-yy)^2+(z2-zz)^2) == c*tdoa12;
    eqn2 = sqrt((x1-xx)^2+(y1-yy)^2+(z1-zz)^2) - sqrt((x3-xx)^2+(y3-yy)^2+(z3-zz)^2) == c*tdoa13;
    eqn3 = sqrt((x1-xx)^2+(y1-yy)^2+(z1-zz)^2) - sqrt((x4-xx)^2+(y4-yy)^2+(z4-zz)^2) == c*tdoa14;

    [x1,y1,z1] = deal(sensor_pos(1,1),sensor_pos(2,1),sensor_pos(3,1));
    [x2,y2,z2] = deal(sensor_pos(1,2),sensor_pos(2,2),sensor_pos(3,2));
    [x3,y3,z3] = deal(sensor_pos(1,3),sensor_pos(2,3),sensor_pos(3,3));
    [x4,y4,z4] = deal(sensor_pos(1,4),sensor_pos(2,4),sensor_pos(3,4));

    sol = solve([eval(eqn1) eval(eqn2) eval(eqn3)],[xx yy zz]);

    solHyp_x1 = matlabFunction(sol.xx(1));
    solHyp_x2 = matlabFunction(sol.xx(2));
    solHyp_y1 = matlabFunction(sol.yy(1));
    solHyp_y2 = matlabFunction(sol.yy(2));
    solHyp_z1 = matlabFunction(sol.zz(1));
    solHyp_z2 = matlabFunction(sol.zz(2));
elseif sum(sensorMask) == 3
    sensor_pos = [rx{sensorInd(1)}.pos rx{sensorInd(2)}.pos rx{sensorInd(3)}.pos];
    
    tdoaPairs(:,2) = sensorInd(2:end); tdoaPairs(:,1) = min(sensorInd);

    syms x1 x2 x3 xx y1 y2 y3 yy c tdoa12 tdoa13 
    eqn1 = sqrt((x1-xx)^2+(y1-yy)^2) - sqrt((x2-xx)^2+(y2-yy)^2) == c*tdoa12;
    eqn2 = sqrt((x1-xx)^2+(y1-yy)^2) - sqrt((x3-xx)^2+(y3-yy)^2) == c*tdoa13;

    [x1,y1] = deal(sensor_pos(1,1),sensor_pos(2,1));
    [x2,y2] = deal(sensor_pos(1,2),sensor_pos(2,2));
    [x3,y3] = deal(sensor_pos(1,3),sensor_pos(2,3));

    sol = solve([eval(eqn1) eval(eqn2)],[xx yy]);

    solHyp_x1 = matlabFunction(sol.xx(1));
    solHyp_x2 = matlabFunction(sol.xx(2));
    solHyp_y1 = matlabFunction(sol.yy(1));
    solHyp_y2 = matlabFunction(sol.yy(2));
else
    sensor_pos = [rx{sensorInd(1)}.pos rx{sensorInd(2)}.pos];
    tdoaPairs = sensorInd;
end

sensorselect = size(sensor_pos,2);
c = physconst('Lightspeed');
% profile off
elapsedTime = toc;
fprintf('Execution time %.3f ms\n', elapsedTime * 1000);
