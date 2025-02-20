%% Emitter Localization using FDoA with Stationary Sensors (Adaptive Grid + Gradient Descent)

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

% This version of the script incorporates Adaptive Grid Refinement and 
% the iterative Quasi-Newton optimization method for FDoA-based emitter 
% localization. The goal is to reduce computational overhead while 
% maintaining accuracy.

% The grid search method is divided into two steps where the search is
% first carried out across the wider span with a lower resolution which
% isolates the region of interest without compromising on the timing and
% computational load. Based on the coarse estimate obtained, a finer grid
% search is initiated in a smaller region. The bounds of this new region
% are defined using the confidence based adaptive region so as not to
% potentially exclude the true position of the emitter from the new search
% area. This method has been derived after experimenting with many other
% fixed bound and non-generalizable techniques.

%% Algorithm for Position Estimation using FDoA (Stage 1)
% 1. Initialize Parameters
% 
% 2. Simulate FDoA measurements with noise between a reference sensor and 
%    the others. 
% 
% 3. Define Search Grid:
%    Create a grid over a 2D search space representing possible emitter 
%    positions.
% 
% 4. Evaluate Cost Function at each coarse grid Point:
%    For each candidate position on the grid:
%    Calculate the range and range rate between the candidate position and 
%    each sensor.
%    Estimate the expected FDoA using the assumed velocity.
%    Compute the cost function as the squared error between the measured 
%    FDoA and estimated FDoA across all CRS pairs.
%    Once the coarse-grid estimate is acquired, a new finer grid is created
%    using the outcome of the coarse-grid search. The finer grid is then
%    searched with a higher resolution to obtain a better estimate of the
%    emitters position
% 
% 5. Implement Unconstrained Nonlinear Optimization using Quasi-Newton:
%    Find the grid point with the minimum cost value.
%    Assign this grid point as the estimated emitter position.

%% Assumptions:
% 1. Single emitter
% 2. A reasonable initial velocity guess
% 3. Sufficient geometric diversity of sensors

%% Caveats
% 1. A poor initial velocity guess can lead to a diverging solution

% Abeer Nasir Chaudhry 
% December 2024



%% Stage 1: Position Estimation with Adaptive Grid Refinement
clc; clear; close all

%% Parameters
c = physconst('LightSpeed');    % Speed of light (m/s)
f0 = 1e9;                       % Carrier frequency (Hz)

sensorPos = [0, 0; % Sensor positions (m)
            10e3, 20e3;
            10e3, -20e3;
            -25e3, 0].';

sensorVel = [0 0;
            0 0;
            0 0;
            0 0].';

Nsensors = size(sensorPos, 2);

emitterPos = [400.35515e3, 250.568e3].'; % True emitter position (m)

v_init = 5e9;                 % Fixed velocity estimate (m/s) for initial approximation
vunitvec_init = [100,100];      % Emitter velocity vector
v_assumed = [120,80];
emitterVel = vunitvec_init.';

%% Simulate FDoA Measurements
fdoaMeas = zeros(1, Nsensors-1);
rangeRate = zeros(1, Nsensors);
range = vecnorm(emitterPos - sensorPos);

for i = 1:Nsensors
    rangeRate(i) = v_init * dot(vunitvec_init, (emitterPos - sensorPos(:,i))) / range(i);
end

fdoaMeas = f0 / c * (rangeRate(2:end) - rangeRate(1));

noise_std = 500;  % Hz

fdoaMeas = fdoaMeas + noise_std * randn(size(fdoaMeas));

%% Adaptive Grid Search
% Initialize Grid Parameters
coarsegridRes = 12.5e3;
[xSpan, ySpan] = deal(500e3, 500e3);

% Check if the source lies outside the search grid
chkSrcPos = abs(emitterPos(1,:)) > xSpan | abs(emitterPos(2,:)) > ySpan;
scidx = find(chkSrcPos==1);
if any(chkSrcPos)
    error('Source %d is outside the search grid. Please expand the search grid or reposition the source',scidx)
end

% profile on

%% Coarse Grid Search
[xCoarse, yCoarse] = meshgrid(-xSpan:coarsegridRes:xSpan,-ySpan:coarsegridRes:ySpan);
coarseGridMap = [xCoarse(:), yCoarse(:)];
coarseCost = zeros(size(coarseGridMap,1), 1);

for idx = 1:length(coarseGridMap)
    pos = coarseGridMap(idx, :).';
    range = vecnorm(pos - sensorPos);
    
    for i = 1:Nsensors

        rangeRate(i) = v_init * dot(v_assumed, (pos - sensorPos(:,i))) / range(i);

    end
    
    fdoaEst = f0 / c * (rangeRate(2:end) - rangeRate(1));
    coarseCost(idx) = sum((fdoaMeas - fdoaEst).^2);
end

[~, idxMin] = min(coarseCost);
coarseBestPos = coarseGridMap(idxMin, :).';

%% Fine Grid Search
% [xFine, yFine] = meshgrid(coarseBestPos(1)-coarsegridRes/2:100:coarseBestPos(1)+coarsegridRes/2, coarseBestPos(2)-coarsegridRes/2:100:coarseBestPos(2)+coarsegridRes/2);

confidenceThreshold = 1.05; 
validIndices = coarseCost <= confidenceThreshold * min(coarseCost);
confidencePoints = coarseGridMap(validIndices, :);

% Extract Fine Grid Bounds
xMin = min(confidencePoints(:,1));
xMax = max(confidencePoints(:,1));
yMin = min(confidencePoints(:,2));
yMax = max(confidencePoints(:,2));

% Validate Bounds
if (xMax - xMin) < 1000
    xMin = coarseBestPos(1) - 5000;
    xMax = coarseBestPos(1) + 5000;
end

if (yMax - yMin) < 1000
    yMin = coarseBestPos(2) - 5000;
    yMax = coarseBestPos(2) + 5000;
end

% Create Fine Grid
[xFine, yFine] = meshgrid(xMin:100:xMax, yMin:100:yMax);

fineGridMap = [xFine(:), yFine(:)];
fineCost = zeros(size(fineGridMap,1), 1);

for idx = 1:length(fineGridMap)
    pos = fineGridMap(idx, :).';
    range = vecnorm(pos - sensorPos);
    
    for i = 1:Nsensors
    
        rangeRate(i) = v_init * dot(v_assumed, (pos - sensorPos(:,i))) / range(i);

    end

fdoaEst = f0 / c * (rangeRate(2:end) - rangeRate(1));
fineCost(idx) = sum((fdoaMeas - fdoaEst).^2);
end

[~, idxMinFine] = min(fineCost);
fineBestPos = fineGridMap(idxMinFine, :).';

%% Gradient Descent Refinement
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off','OptimalityTolerance',1e-9);
optimizedPos = fminunc(@(pos) costFunction(pos, sensorPos, fdoaMeas, f0, c, v_init, vunitvec_init),fineBestPos,options);

fprintf(['True Position: (%.3f, %.3f) km\n' ...
         'Coarse Estimate: (%.3f, %.3f) km\n'...
         'Fine Estimate: (%.3f, %.3f) km\n'...
         'Optimized Estimate: (%.3f, %.3f) km\n'], emitterPos(1)/1e3,emitterPos(2)/1e3, ...
                                                   coarseBestPos(1)/1e3, coarseBestPos(2)/1e3, ...
                                                   fineBestPos(1)/1e3, fineBestPos(2)/1e3, ...
                                                   optimizedPos(1)/1e3, optimizedPos(2)/1e3);
% profile viewer
%% Visualization
hold on; grid on;
scatter(sensorPos(1,:),sensorPos(2,:),'MarkerFaceColor','b','MarkerEdgeColor','k','Marker','o','DisplayName','Sensor Network');
scatter(emitterPos(1,:),emitterPos(2,:),'MarkerFaceColor','r','MarkerEdgeColor','k','Marker','^','DisplayName','True Emitter Location');
scatter(coarseBestPos(1,:),coarseBestPos(2,:),'MarkerFaceColor','y','MarkerEdgeColor','k','Marker','o','DisplayName','Coarse Position Esimate');
scatter(optimizedPos(1,:),optimizedPos(2,:),'MarkerFaceColor','g','MarkerEdgeColor','k','Marker','diamond','DisplayName','Optimized Position Estimate');
plot([coarseBestPos(1) optimizedPos(1)],[coarseBestPos(2) optimizedPos(2)],'LineStyle',':','DisplayName','Optimization')

% Plot the velocity vector
quiver(emitterPos(1,:),emitterPos(2,:),vunitvec_init(1),vunitvec_init(2),4,'MaxHeadSize',2,'LineWidth',1.5,'DisplayName','Velocity Vector','Color','r')

for i = 2:(size(sensorPos,2))
    label = ('  S'+ string(i));
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
    vdiff = CalcDopDiff(emitterPos,emitterVel,sensorPos(:,1),sensorVel(:,1),sensorPos(:,i),sensorVel(:,i),c);
    xy_isodop = drawIsodop(sensorPos(:,1),sensorVel(:,1),sensorPos(:,i),sensorVel(:,i),emitterVel,vdiff,1000,max(max(emitterPos))*1.5);
    hh = plot(xy_isodop(1,:),xy_isodop(2,:),'DisplayName','Line of Constant FDOA','LineStyle','--','LineWidth',0.1,'Color','k');
    excludeFromLegend(hh);

end

xlim([-xSpan xSpan]);ylim([-ySpan ySpan])
legend;

legend;
title('Adaptive Grid and  QN optimization for FDoA Localization');
xlabel('Cross-Range (m)');
ylabel('Down-Range (m)');

%% Error Analysis
posError = norm(emitterPos - optimizedPos);
fprintf('Position Error: %.2f m\n', posError);


%% Cost Function Definition
function cost = costFunction(pos, sensorPos, fdoaMeas, f0, c, v_init, vunitvec_init)
    range = vecnorm(pos - sensorPos);
    for i = 1:size(sensorPos, 2)
        rangeRate(i) = v_init * dot(vunitvec_init, (pos - sensorPos(:,i))) / range(i);
    end
    fdoaEst = f0 / c * (rangeRate(2:end) - rangeRate(1));
    cost = sum((fdoaMeas - fdoaEst).^2);
end

