%% Platfomr positions
clear;clc
sensor_pos = [0.1e3 0.25e3 0.5e3;
              0.05e3 0.25e3 0e3];
nSensors = size(sensor_pos,2);

 source = [10e3;
               90e3];

% Central sensor pairing (ref sensor = s_n)
tdoaPairs(:,2) = (2:(nSensors))';
tdoaPairs(:,1) = 1;


ranges = vecnorm(sensor_pos - source);
% theta = 90;
for i = 1:size(tdoaPairs,1)
        theta = atan2d((sensor_pos(2,i+1) - sensor_pos(2,1)),((sensor_pos(1,i+1) - sensor_pos(1,1))))+90;
        iso{i} = isochrone_func(sensor_pos(:,1),sensor_pos(:,i+1),(ranges(1)-ranges(i+1)), ...
            100e3,5e3,theta);
end

figure;
for i = 1:size(tdoaPairs,1)
    plot(iso{i}(1,:),iso{i}(2,:))
% hold on;plot(iso{2}(1,:),iso{2}(2,:))
    hold on
    scatter(sensor_pos(1,:),sensor_pos(2,:))
    scatter(source(1),source(2))
end

% plot(-iso{2}(1,:),-iso{2}(2,:)); plot(-iso{1}(1,:),-iso{1}(2,:)); plot(-iso{3}(1,:),-iso{3}(2,:))
