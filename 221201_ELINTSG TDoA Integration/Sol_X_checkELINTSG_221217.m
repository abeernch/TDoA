%% TDOA Analytical Solution Script for 3D position estimation for non-cooperative emitter using Hyperboloids

% It is a 3D problem with 3 unknowns in the presence of 'N' sensors, where 
% N > 3.

% A centralized reference sensor (CRS) approach is used to solve the
% equations of the hyperboloids with S1 as the reference sensor. The
% drawback of ths CRS approach is the inability to robustly distinguish 
% non-unique (GHOST) targets by having multiple solutions using different
% baselines. 

% The advantage is limited computation requirements. This is an alternative
% approach to decentralized reference sensor (DRS) network.

% This is an extension of the test/demo script
% "analyticalSol_3D_Hyperboloid_221118". This script aims to benchmark the
% computational complexity of the solution process for 'M' TOA readings.
% The solution has been made efficient by solving the analytical equation
% only once partially based on the available data, i.e., sensor positions
% and 'c'. The partially solved equations for the variables x, y and z are
% saved as function handles with t1, t2, t3 and t4 as inputs to the
% handles.
% The solutions for 500 TOA readings are obtained in under 100ms


% Assumptions: 
% 1. No noise 2. No multipath effects


% Author: Abeer Nasir Chaudhry®
% 221124

% clear; clc

graph = 0;

%% Initialize variables and CAS equations
syms x1 x2 x3 x4 x y1 y2 y3 y4 y z1 z2 z3 z4 z t1 t2 t3 t4 t c
eqn1 = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2) - sqrt((x2-x)^2 + (y2-y)^2 + (z2-z)^2)  == c*(t1-t2);
eqn2 = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2) - sqrt((x3-x)^2 + (y3-y)^2 + (z3-z)^2)  == c*(t1-t3);
eqn3 = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2) - sqrt((x4-x)^2 + (y4-y)^2 + (z4-z)^2)  == c*(t1-t4);


%% Genrerate TOAs given scenario
sensor_pos = [rx{1}.pos rx{2}.pos rx{3}.pos rx{4}.pos];

% Duration of emission
% t = linspace(0,360,10000);

% Constant emitter platform velocity
% vtx = [50 50 50];
% r = 100e3;
for i = 1:length(tx{1}.trajPos.')
    em_pos = tx{1}.trajPos.';
    emitter_pos(1:3,i) = em_pos(1:3,i);

    [nSensor, ~, ~, ~, ~, reported_toa_int, range] = mod_ToaGenerator(sensor_pos,emitter_pos(1:3,i),graph);
    toa_arr(:,i) = reported_toa_int+ (rand(4,1)*(5e-9))-(0e-9); 
end

%% Assign parameters to equation variables
c = 3e8;
[x1,y1,z1]       = deal(sensor_pos(1,1),sensor_pos(2,1),sensor_pos(3,1));
[x2,y2,z2]       = deal(sensor_pos(1,2),sensor_pos(2,2),sensor_pos(3,2));
[x3,y3,z3]       = deal(sensor_pos(1,3),sensor_pos(2,3),sensor_pos(3,3));
[x4,y4,z4]       = deal(sensor_pos(1,4),sensor_pos(2,4),sensor_pos(3,4));

% Solve partial equation and create inline functions
sol = solve([eval(eqn1) eval(eqn2) eval(eqn3)],[x y z],"IgnoreAnalyticConstraints",true);
solHyp_x1 = matlabFunction(sol.x(1));
solHyp_x2 = matlabFunction(sol.x(2));
solHyp_y1 = matlabFunction(sol.y(1));
solHyp_y2 = matlabFunction(sol.y(2));
solHyp_z1 = matlabFunction(sol.z(1));
solHyp_z2 = matlabFunction(sol.z(2));

%% Intialize for loop 
tic

for i = 1:length(toa_arr)
    
    [t1,t2,t3,t4] = deal(toa_arr(1,i),toa_arr(2,i),toa_arr(3,i),toa_arr(4,i));
    
    %% Evaluate equations
    sol_x(1:2,i) = [solHyp_x1(t1,t2,t3,t4); solHyp_x2(t1,t2,t3,t4)];
    sol_y(1:2,i) = [solHyp_y1(t1,t2,t3,t4); solHyp_y2(t1,t2,t3,t4)];
    sol_z(1:2,i) = [solHyp_z1(t1,t2,t3,t4); solHyp_z2(t1,t2,t3,t4)];
    disp(i)
end

toc

figure;p = plot3(em_pos(1,:),em_pos(2,:),em_pos(3,:),'LineWidth',3,'DisplayName','Actual Platform Trajectory');
hold on
sols = scatter3(sol_x(1,:),sol_y(1,:),sol_z(1,:),'red','o','LineWidth',0.1,'DisplayName','Measured Trajectory','MarkerEdgeAlpha',0.25);
sols2 = scatter3(sol_x(2,:),sol_y(2,:),sol_z(2,:),'green','o','LineWidth',0.1,'DisplayName','Measured Trajectory','MarkerEdgeAlpha',0.25);
s_pos =  scatter3(sensor_pos(1,:),sensor_pos(2,:),sensor_pos(3,:),'blue','filled','o','LineWidth',0.1,'DisplayName','Sensor Network','MarkerEdgeAlpha',1);
legend(p.DisplayName,sols.DisplayName,sols2.DisplayName,s_pos.DisplayName,'Location','northeast')
grid on
hold off
figure
p = plot3(em_pos(1,:),em_pos(2,:),em_pos(3,:),'LineWidth',3,'DisplayName','Actual Platform Trajectory');
hold on
sols = scatter3(sol_x(1,:),sol_y(1,:),sol_z(1,:),'red','o','LineWidth',0.1,'DisplayName','Solution 1','MarkerEdgeAlpha',0.25);
sols2 = scatter3(sol_x(2,:),sol_y(2,:),sol_z(2,:),'green','o','LineWidth',0.1,'DisplayName','Solution 2','MarkerEdgeAlpha',0.25);
legend(p.DisplayName,sols.DisplayName,sols2.DisplayName,'Location','northeast')

%% PLOTTING SOLUTIONS
%% SOLUTION 1
sol1 = [sol_x(1,:); sol_y(1,:); sol_z(1,:)];
err_x1 = abs(emitter_pos(1,:) - sol1(1,:));
err_y1 = abs(emitter_pos(2,:) - sol1(2,:));
err_z1 = abs(emitter_pos(3,:) - sol1(3,:));

figure; 
plot(1:length(err_x1),err_x1, ...
    1:length(err_y1),err_y1, ...
    1:length(err_z1),err_z1, ...
    1:length(err_z1),vecnorm(emitter_pos - sol1), ...
    'LineWidth',1.5); title('Error in Solution "1"'); xlabel('TOA Sample #'); ylabel('Error (meters)')
legend('Error in X','Error in Y','Error in Z', 'Absolute Range Error')

figure; 
subplot(4,1,1); scatter(1:length(err_z1),vecnorm(emitter_pos - sol1),15,'filled'); title('Absolute Range Error'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(4,1,2); scatter(1:length(err_x1),err_x1,15,'filled'); title('Error in X'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(4,1,3); scatter(1:length(err_y1),err_y1,15,'filled'); title('Error in Y'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(4,1,4); scatter(1:length(err_z1),err_z1,15,'filled'); title('Error in Z'); xlabel('TOA Sample #'); ylabel('Error (meters)')

%% SOLUTION 2
sol2 = [sol_x(2,:); sol_y(2,:); sol_z(2,:)];
err_x2 = (em_pos(1,:) - sol2(1,:));
err_y2 = (em_pos(2,:) - sol2(2,:));
err_z2 = (em_pos(3,:) - sol2(3,:));

figure; 
plot(1:length(err_x1),err_x2, ...
    1:length(err_y1),err_y2, ...
    1:length(err_z1),err_z2, ...
    1:length(err_z1),vecnorm(em_pos - sol2), ...
    'LineWidth',1.5); title('Error in Solution "2"'); xlabel('TOA Sample #'); ylabel('Error (meters)')
legend('Error in X','Error in Y','Error in Z', 'Absolute Range Error')

figure; 
subplot(4,1,1); plot(1:length(err_z1),vecnorm(em_pos - sol2)); title('Absolute Range Error'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(4,1,2); plot(1:length(err_x1),abs(err_x1)); title('Error in X'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(4,1,3); plot(1:length(err_y1),abs(err_y1)); title('Error in Y'); xlabel('TOA Sample #'); ylabel('Error (meters)')
subplot(4,1,4); plot(1:length(err_z1),abs(err_z1),LineWidth=0.1); title('Error in Z'); xlabel('TOA Sample #'); ylabel('Error (meters)')
%% Compute RMSE for both solutions

% Solution 1
x1_rmse = sqrt(sum(abs(err_x1.^2))/length(err_x1))
y1_rmse = sqrt(sum(abs(err_y1.^2))/length(err_y1))
z1_rmse = sqrt(sum(abs(err_z1.^2))/length(err_z1))

% Solution 2
x2_rmse = sqrt(sum(abs(err_x2.^2))/length(err_x2))
y2_rmse = sqrt(sum(abs(err_y2.^2))/length(err_y2))
z2_rmse = sqrt(sum(abs(err_z2.^2))/length(err_z2))