%% TDOA Analytical Solution Script for 3D position estimation for non-cooperative emitter using Hyperboloids

% It is a 3D problem with 3 unknowns in the presence of 'N' sensors, where 
% N > 3.

% A centralized reference sensor (CRS) approach is used to solve the
% equations of the hyperboloids with S1 as the reference sensor. The
% drawback of ths CRS approach is the inability to robustly distinguish 
% non-unique (GHOST) targets. 

% The advantage is limited computation requirements. This is an alternative
% approach to decentralized reference sensor (DRS) network.

% This is an extension of the test/demo script
% "analyticalSol_3D_Hyperboloid_221118". This script aims to benchmark the
% computational complexity of the solution process for 'M' TOA readings.


% Assumptions: 
% 1. No noise 2. No multipath effects


% Author: Abeer Nasir Chaudhry®
% 221124

clear;clc

graph = 0;

%% Initialize variables and CAS equations
syms x1 x2 x3 x4 x y1 y2 y3 y4 y z1 z2 z3 z4 z t1 t2 t3 t4 t c
eqn1 = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2) - sqrt((x2-x)^2 + (y2-y)^2 + (z2-z)^2)  == c*(t1-t2);
eqn2 = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2) - sqrt((x3-x)^2 + (y3-y)^2 + (z3-z)^2)  == c*(t1-t3);
eqn3 = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2) - sqrt((x4-x)^2 + (y4-y)^2 + (z4-z)^2)  == c*(t1-t4);


%% Genrerate TOAs given scenario
sensor_pos = [0 -10e3 10e3 00e3; 
              0  10e3 10e3 -15e3;
              0  00e3 00e3  00e3];

% Duration of emission
t = linspace(0.1,100,50);

% Constant emitter platform velocity
vtx = [10 50 2];

for i = 1:length(t)
    emitter_pos = [t(i)*vtx(1); t(i)*vtx(2);t(i)*vtx(3)];
    [nSensor, ~, ~, ~, ~, reported_toa, range] = mod_ToaGenerator(sensor_pos,emitter_pos,graph);
    toa_arr(:,i) = reported_toa; 
end

%% Assign parameters to equation variables
c = 3e8;
[x1,y1,z1]       = deal(sensor_pos(1,1),sensor_pos(2,1),sensor_pos(3,1));
[x2,y2,z2]       = deal(sensor_pos(1,2),sensor_pos(2,2),sensor_pos(3,2));
[x3,y3,z3]       = deal(sensor_pos(1,3),sensor_pos(2,3),sensor_pos(3,3));
[x4,y4,z4]       = deal(sensor_pos(1,4),sensor_pos(2,4),sensor_pos(3,4));
%% Intialize for loop 
tic

for i = 1:length(toa_arr)
    [t1,t2,t3,t4] = deal(toa_arr(1,i),toa_arr(2,i),toa_arr(3,i),toa_arr(4,i));
    disp([i])
    %% Evaluate equations
    s1 = eval(eqn1);
    s2 = eval(eqn2);
    s3 = eval(eqn3);
    [xsol, ysol, zsol] = solve([s1 s2 s3],[x y z]);
    %% Extract soln
    x_est(1:2,i) = (eval(xsol));
    y_est(1:2,i) = (eval(ysol));
    z_est(1:2,i) = (eval(zsol)); 
end

toc