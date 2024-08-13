% Initialize Detector/Source Locations
x_sensor1 = [pos_sensor(1,1); pos_sensor(2,1)];
x_sensor2 = [pos_sensor(1,2);pos_sensor(2,2) ];
x_sensor3 = [pos_sensor(1,3);pos_sensor(2,3)];
x_sensor4 = [pos_sensor(1,4);pos_sensor(2,4)];
x_source = [pos_emitter(1,1); pos_emitter(2,1)];
% Compute Ranges
r1 = utils.rng(x_sensor1,x_source);
r2 = utils.rng(x_sensor2,x_source);
r3 = utils.rng(x_sensor3,x_source);
r4 = utils.rng(x_sensor4,x_source);
% Find Isochrones
xiso1 = tdoa.drawIsochrone(x_sensor1,x_sensor4,r4-r1,10e3,50e3);
xiso2 = tdoa.drawIsochrone(x_sensor2,x_sensor4,r4-r2,10e3,100e3);
xiso3 = tdoa.drawIsochrone(x_sensor3,x_sensor4,r4-r3,10e3,100e3);


% Draw Figure
fig1b = figure();hold on;
% Isochrones
plot(xiso1(1,:),xiso1(2,:),'k:','DisplayName','Isochrone');
hiso2=plot(xiso2(1,:),xiso2(2,:),'k:');
utils.excludeFromLegend(hiso2);
hiso3=plot(xiso3(1,:),xiso3(2,:),'k:');
utils.excludeFromLegend(hiso3);
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