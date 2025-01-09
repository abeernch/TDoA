sensorPos = [00e3	-20e3	10e3	10e3
             0	  0 	   20e3	 -20e3];
refSen = 4;
timingErr = 1e-9;
res = 2e3;
gridSpan = [-400e3,400e3; -400e3, 400e3];
[X,Y,ggg] = gdop2Dsd(sensorPos,timingErr, res,gridSpan,refSen);

contourLevelsxy = 100:100:2000;
figure;
scatter(sensorPos(1,:),sensorPos(2,:))
hold on
contour(X,Y,ggg,contourLevelsxy,'LineStyle','--','LineColor','k')