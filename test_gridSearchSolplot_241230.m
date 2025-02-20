%% Emitter Localization using FDoA with Stationary Sensors

% This is an approach that is under consideration for the FDoA-based
% localization of a scenario including a non-stationary emitter with
% stationary sensors. 
% This is a two-step solver. The process is divided into two stages where
% the estimation of the position of the emitter is first calculated with an
% intial guess for the emitters velocity. With the variables reduced down
% to 3 (x,y,z) instead of 6 (x,y,z,vx,vy,vz), a CRS-based Solution of the
% 3 baseline FDoA equations yield the solution which is the location of the
% emitter. Once the position has been estimated, the solution can be used
% again by susbstituting the position to determine the velocity of the 
% emitter.

% This script incorparates the Grid Search solution and the plotting of the
% contours.

%% Algorithm for Position Estimation using FDoA (Stage 1)
% 1. Initialize Parameters:
% 2. Simulate FDoA measurements with noise between a reference sensor and the others.
% 
% 3. Define Search Grid:
%   Create a grid over a 2D search space representing possible emitter positions.
% 
% 4. Evaluate Cost Function at Each Grid Point:
%   For each candidate position on the grid:
%   Calculate the range and range rate between the candidate position and each sensor.
%   Estimate the expected FDoA using the assumed velocity.
%   Compute the cost function as the squared error between the measured FDoA and estimated FDoA across all sensor pairs.
% 
% 5. Identify Minimum Cost Point:
%   Find the grid point with the minimum cost value.
%   Assign this grid point as the estimated emitter position.

%% Assumptions:
% 1. Single emitter
% 2. A reasonable initial velocity guess
% 3. Sufficient geometric diversity of sensors

%% Caveats
% 1. A poor initial guess can lead to a diverging solution

% Abeer Nasir Chaudhry 
% December 2024

%% Stage 1: Position Estimation with Fixed Velocity
clc; clear; close all

%% Parameters
c = physconst('LightSpeed');    % Speed of light (m/s)
f0 = 1e9;                       % Carrier frequency (Hz)

sensorPos = [00, 00; % Sensor positions (m)
            10e3, 20e3;
            10e3, -20e3;
           -10e3, 0].';

sensorVel = [0 0; 
             0 0; 
             0 0; 
             0 0].';

Nsensors = size(sensorPos, 2);

emitterPos = [264.27e3, -95.1094e3].'; % True emitter position (m)

v_init = 4000;                   % Fixed velocity estimate (m/s) for initial approximation
vunitvec_init = [1,1];          % Emitter traveling along the y-axis
emitterVel =vunitvec_init.';
profile on

%% Simulate FDoA Measurements
fdoaMeas = zeros(1, Nsensors-1);                                                        % Pre-allocation     
rangeRate = zeros(1,Nsensors);

range = vecnorm(emitterPos - sensorPos);                                                % Calculate range from each sensor to the emitter

for i = 1:Nsensors
    rangeRate(i) = v_init*dot(vunitvec_init, (emitterPos - sensorPos(:,i))) / range(i); % Calculate the Projected Doppler on each sensor  
end

fdoaMeas = f0 / c * (rangeRate(2:end) - rangeRate(1));                                  % Compute the FDoA


noise_std = 2; % Add measurement noise
fdoaMeas = fdoaMeas + noise_std * randn(size(fdoaMeas));

%% Position Estimation via Iso-Doppler Curves
% A search grid is created over a plausible range of x and y coordinates.
% At each grid point, the FDoA values are predicted assuming a constant velocity direction 
% A cost function calculates the difference between the simulated noisy FDoA values and the predicted FDoA values at each grid point

% Initialize grid and grid parameters
gridRes = 12.5e3;
[xSpan, ySpan] = deal(500e3, 500e3);

[xGrid, yGrid] = meshgrid(-xSpan:gridRes:xSpan, -ySpan:gridRes:ySpan);

gridMap = [xGrid(:),yGrid(:)];

% Check if the source lies outside the search grid
chkSrcPos = abs(emitterPos(1,:)) > abs(min(gridMap(1))) | abs(emitterPos(2,:)) > abs(min(gridMap(2)));
scidx = find(chkSrcPos==1);
if any(chkSrcPos)
    error('Source %d is outside the search grid. Please expand the search grid or reposition the source',scidx)
end

% Initialize cost variable
cost = zeros(size(xGrid));

% Begin the grid search and compute the FDoA for each possible x,y. This is
% then compared with the actual FDoA to generate the cost at that point.
% The loop is carried out for each point in the grid
% (size(xGrid,1)*size(yGrid,1))

for idx = 1:length(gridMap)
    pos = gridMap(idx, :).'; % Candidate emitter position
    costVal = 0;

    range = vecnorm(pos - sensorPos);

    for i = 1:Nsensors
        
        rangeRate(i) = v_init * dot(vunitvec_init, (pos - sensorPos(:,i))) / range(i); 
        
    end

    fdoaEst = f0 / c * (rangeRate(2:end) - rangeRate(1));
    costVal = sum((fdoaMeas - fdoaEst).^2);     % Evaluate the cost from FDoA estimated from each grid position against the measured FDoA
    cost(idx) = costVal;
end

% Find the grid point with minimum cost
[minCostval,minCostidx] = min(cost(:));

EstPos = gridMap(minCostidx,:).';
fprintf('Estimated Position: (%.2f, %.2f) m\n', EstPos(1), EstPos(2));
profile viewer

%% Plot
%% Visualization
% figure;
hold on; grid on;
scatter(sensorPos(1,:),sensorPos(2,:),'MarkerFaceColor','b','MarkerEdgeColor','k','Marker','o','DisplayName','Sensor Network'); % Plot sensor positions
text(sensorPos(1,1)-.2,sensorPos(2,1)-1.2,'S1','FontSize',10);
scatter(emitterPos(1,:),emitterPos(2,:),'MarkerFaceColor','r','MarkerEdgeColor','k','Marker','^','DisplayName','True Emitter Location');
scatter(EstPos(1,:),EstPos(2,:),'MarkerFaceColor','g','MarkerEdgeColor','k','Marker','diamond','DisplayName','Estimated Emitter Position');

% Plot the velocity vector
quiver(emitterPos(1,:),emitterPos(2,:),vunitvec_init(1),vunitvec_init(2),4,'MaxHeadSize',2,'LineWidth',1.5,'DisplayName','Velocity Vector','Color','r')

for i = 2:(size(sensorPos,2))
    label = ('  S'+ string(i));
    % scatter(sensorPos(1,:),sensorPos(2,:),'DisplayName','Sensors','Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k');
    text(sensorPos(1,i)-.2,sensorPos(2,i)-.2,label,'FontSize',10);

    % Draw velocity arrows for emitter
    for j = 1:size(emitterPos,2)
      % srcmtn = quiver(emitterPos(1,j),emitterPos(2,j),(emitterVel(1,j)/4), (emitterVel(2,j)/4),'MaxHeadSize',4,'Color','k');
      % excludeFromLegend(srcmtn);
    end

    % Draw velocity arrows for sensors
   % senmtn = quiver(sensorPos(1,i),sensorPos(2,i),(sensorVel(1,i)/4), (sensorVel(2,i)/4),'MaxHeadSize',4, 'Color','k','LineStyle','-');
   % excludeFromLegend(senmtn);

    % Draw isodoppler line S12
    vdiff = CalcDopDiff(emitterPos,emitterVel,sensorPos(:,1),sensorVel(:,1),sensorPos(:,i),sensorVel(:,i),3e8);
    xy_isodop = drawIsodop(sensorPos(:,1),sensorVel(:,1),sensorPos(:,i),sensorVel(:,i),emitterVel,vdiff,1000,max(max(emitterPos))*1.5);
    hh = plot(xy_isodop(1,:),xy_isodop(2,:),'DisplayName','Line of Constant FDOA','LineStyle','--','LineWidth',0.1,'Color','k');
end

excludeFromLegend(hh);
xlim([-xSpan xSpan]);ylim([-ySpan ySpan])
legend;
title('Two Step Emitter Position Estimation via FDoA');
xlabel('Cross-Range (km)');
ylabel('Down-Range (km)');

%% Error Analysis
posError = norm(emitterPos - EstPos);
fprintf('Position Error:%0.2f ', posError);

fprintf('Localization complete.\n');

%% Script Log
% 1. (Dec 22, 2024). File created.

% 2. The velocity estimate that is used to decouple position and velocity
% in the 2-step Grid Search Method allows a good enough baseline for the
% Doppler projections from the source to the sensors.
% The iso-Doppler contours are geometry driven with a focus on
% emitter-sensor relative geometry. It can also be observed from the
% equation for the range-rate difference that the velocity acts as a
% scaling factor only. Unless the velocity estimate is off by several
% orders, or the source is located outside the search grid, the algorithm
% works.