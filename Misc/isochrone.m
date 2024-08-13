% Initialize Detector/Source Locations
clear;clc;close all
grid on
axis([-100 100 -100 100])
[x,y] = ginput(3);
close;
x_sensor1 = [x(1); y(1)];
x_sensor2 = [x(2); y(2)];
x_sensor3 = [x(3); y(3)];

grid on
axis([-100 100 -100 100])
[x,y] = ginput(1);
close;
x_source = [x; y];
% Compute Ranges
r1 = utils.rng(x_sensor1,x_source);
r2 = utils.rng(x_sensor2,x_source);
r3 = utils.rng(x_sensor3,x_source);
% Find Isochrones
xiso1 = tdoa.drawIsochrone(x_sensor1,x_sensor2,r2-r1,200,100);
xiso2 = tdoa.drawIsochrone(x_sensor2,x_sensor3,r3-r2,200,100);
xiso3 = tdoa.drawIsochrone(x_sensor1,x_sensor3,r3-r1,200,100);


% Draw Figure
fig1b = figure();hold on;
% Isochrones
plot(xiso1(1,:),xiso1(2,:),'k:','DisplayName','Isochrone');
hiso2=plot(xiso2(1,:),xiso2(2,:),'k:');
hiso3=plot(xiso3(1,:),xiso3(2,:),'k:');

utils.excludeFromLegend(hiso2);
% Isochrone Labels
text(mean([x_sensor1(1),x_sensor2(1)]),mean([x_sensor1(2),x_sensor2(2)])-.2,'TDOA_{1,2}');
text(mean([x_sensor2(1),x_sensor3(1)])+.3,mean([x_sensor2(2),x_sensor3(2)]),'TDOA_{2,3}');
% Position Markers
hiso2=plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'k-','LineWidth',1);
utils.excludeFromLegend(hiso2);
plot([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'ko','DisplayName','Sensors');
plot(x_source(1),x_source(2),'k^','MarkerSize',8,'DisplayName','Transmitter');
% Position Labels
text(x_sensor1(1)+.05,x_sensor1(2)-.1,'S_1');
text(x_sensor2(1)+.05,x_sensor2(2)-.1,'S_2');
text(x_sensor3(1)+.05,x_sensor3(2)-.1,'S_3');
% text(x_source(1)+.05,x_source(2)+.05,'$T_1$');
% Adjust Axes
xlim([min([xiso1(1,:) xiso2(1,:)])-15 max([xiso1(1,:) xiso2(1,:)])]);
ylim([min([xiso1(2,:) xiso2(2,:)])-15 max([xiso1(2,:) xiso2(2,:)])]);
legend('Location','SouthEast');
% utils.setPlotStyle(gca,{'clean','widescreen','tight'});
% utils.exportPlot(fig1b,[prefix '1b']);
