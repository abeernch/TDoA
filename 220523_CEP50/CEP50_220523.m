% clear;close all;clc

% Define Sensor Positions
figure
axis([-100 100 -100 100])
grid on
[x,y] = ginput();
close
nSensors = length(x);
sensor_pos = zeros(2,length(x));
sensor_pos(1,:) = x*1e3;
sensor_pos(2,:) = y*1e3;
% Define Sensor Performance
timingError = 30e-9;
% Create the error covariance matrix (time of arrival)
Ctoa = timingError^2*eye(nSensors);

% Define source positions
grd_size = 1e2; % Grid resolution
x_ax = linspace(-200e3,200e3,grd_size);
y_ax = linspace(-200e3,200e3,grd_size);
[xx,yy] = ndgrid(x_ax,y_ax);
source_posgrid = [xx(:) yy(:)]';

% Compute CRLB and CEP50
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = computeCRLB(sensor_pos,source_posgrid,Ctoa); % Ndim x Ndim x M^2
cep50 = reshape(computeCEP50(crlb),[grd_size,grd_size]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning

% Set up contours
contourLevels=[.1,1,2,3,5,10,15,20];
contourLevelsLabel = [.1,1,2,3,5,10,15,20];

%% Draw Figure
fig = figure();hold on;

plot(sensor_pos(1,:)/1e3,sensor_pos(2,:)/1e3,'o','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);
excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;

% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
% setPlotStyle(gca,{'equal','tight'});
