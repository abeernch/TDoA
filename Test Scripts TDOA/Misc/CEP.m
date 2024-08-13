% Define Sensor Positions
baseline = 10e3;
nSensors = 10;
thSensors = linspace(0,2*pi,nSensors+1) +pi/2; % add an extra sample, will be ignored
x_sensor = baseline* [ cos(thSensors(1:end-1));sin(thSensors(1:end-1))];
% Define Sensor Performance
timingError = 1e-7;
Ctoa = timingError^2*eye(nSensors); % utilities now resample cov matrix with ref_idx
% Ctdoa = timingError^2 * (1 + eye(nSensors-1));
% Define source positions
M = 501;
xvec = linspace(-100,100,M)*1e3;
yvec = linspace(-100,100,M)*1e3;
[xx,yy] = ndgrid(xvec,yvec);
x_source = [xx(:) yy(:)]';
% Compute CRLB
warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
crlb = tdoa.computeCRLB(x_sensor,x_source,Ctoa); % Ndim x Ndim x M^2
cep50 = reshape(utils.computeCEP50(crlb),[M,M]);
warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning
% Set up contours
contourLevels=[.1,1,5,10,50,100,1000];
contourLevelsLabel = [.1,1,5,10,50,100];
% Draw Figure
fig6a = figure();hold on;
%ax=subplot(2,1,1)
plot(x_sensor(1,:)/1e3,x_sensor(2,:)/1e3,'o','DisplayName','Sensors','LineWidth',1);
[cp,hiso2]=contour(xx/1e3,yy/1e3,cep50/1e3,contourLevels,'LineColor','k');
clabel(cp,hiso2);%contourLevelsLabel,fmt='%1.0f',fontsize=10,rightside_up=True)
utils.excludeFromLegend(hiso2);
legend('Location','NorthEast');
grid off;
% Adjust the Display
xlabel('Cross-range [km]');ylabel('Down-range [km]');
% utils.setPlotStyle(gca,{'equal','tight'});
% utils.exportPlot(fig6a,[prefix '6a']);