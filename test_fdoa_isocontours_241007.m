%% 2D FDoA iso doppler contours

% The scripts uses the emitter and sensor network positions and velocities
% to calculate the projection of the emitter(s) velocity on to each sensor followed
% by the difference of the experienced doppler shift at each sensor with
% respect to the central reference sensor (CRS) (idx = 1). The differences in doppler
% shifts w.r.t to the CRS is then used to compute and plot the isocontours
% that end up intersecting at the position of the emitting source.

clear;clc
%% iso Doppler Contours

emitterPos = [200 100].' ; % Source Position 

emitterVel = [1 1].'; % meters per sec

sensorPos = [00, 00; % Sensor positions (m)
            10, 10;
            10, -10;
           -10, 0].';
sensorVel = [0 0; 
             0 0; 
             0 0; 
             0 0].';
% figure;
% Draw Geometry
scatter(emitterPos(1,:),emitterPos(2,:),60,'DisplayName','Transmitter','Marker','^','MarkerFaceColor','r','MarkerEdgeColor','k');
hold on
label = [];

% Initializing graphics for reference sensor
ax = gca;
ax.Units ='normalized';
text(sensorPos(1,1)-.2,sensorPos(2,1)-.2,'S1','FontSize',10);
sen1mtn = quiver(sensorPos(1,1),sensorPos(2,1),(sensorVel(1,1)/4), (sensorVel(2,1)/4),'MaxHeadSize',4, 'Color','k');
excludeFromLegend(sen1mtn)

for i = 2:(size(sensorPos,2))
    label = ('  S'+ string(i));
    scatter(sensorPos(1,:),sensorPos(2,:),'DisplayName','Sensors','Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k');
    text(sensorPos(1,i)-.2,sensorPos(2,i)-.2,label,'FontSize',10);

    % Draw velocity arrows for emitter
    for j = 1:size(emitterPos,2)
      srcmtn = quiver(emitterPos(1,j),emitterPos(2,j),(emitterVel(1,j)/4), (emitterVel(2,j)/4),'MaxHeadSize',4,'Color','k');
      excludeFromLegend(srcmtn);
    end

    % Draw velocity arrows for sensors
   senmtn = quiver(sensorPos(1,i),sensorPos(2,i),(sensorVel(1,i)/4), (sensorVel(2,i)/4),'MaxHeadSize',4, 'Color','k','LineStyle','-');
   excludeFromLegend(senmtn);

    % Draw isodoppler line S12
    vdiff = CalcDopDiff(emitterPos,emitterVel,sensorPos(:,1),sensorVel(:,1),sensorPos(:,i),sensorVel(:,i),3e8);
    xy_isodop = drawIsodop(sensorPos(:,1),sensorVel(:,1),sensorPos(:,i),sensorVel(:,i),emitterVel,vdiff,1000,max(max(emitterPos))*1.5);
    hh = plot(xy_isodop(1,:),xy_isodop(2,:),'DisplayName','Line of Constant FDOA','LineStyle','--','LineWidth',0.1,'Color','k');
    
end 

xlabel('Cross-Range (km)'); ylabel('Down-Range (km)'); title('FDoA Isodoppler Contours '+ string(size(sensorPos,2)) + ' Sensors')
legend('Source', 'Sensors', 'FDoA Solution')
excludeFromLegend(hh);
