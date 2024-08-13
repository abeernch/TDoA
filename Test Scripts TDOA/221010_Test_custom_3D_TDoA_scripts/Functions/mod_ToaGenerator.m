function [nSensors, nEmitter, pos_sensor, pos_emitter, emitter_init, reported_toa, range] = mod_ToaGenerator(sensor_pos,emitter_pos,graph)
%% Generates TOA at sensors of a network for multilateration testing.
% - The positions for the sensors and the emitter can be defined
% interactively within the function, or they can be manually input if the
% interactive plotting is disabled by setting the 'graph' input as false
% - Other inputs are total time and PRI. These are used to create an array
% of TOAs based on a preselected PRI for all sensors for the entered total
% time
% - Based on the interactively (or otherwise) placed platforms, the ranges
% are calculated and converted to TOA at each sensor. 
% - The is adapated to cater for any no of sensors and emitters
% - At the output, the platform positions, toa array and no of platforms
%  are acquired.

% Author: Abeer Nasir Chaudhry®

if graph == 1
    %% Initialize Sensor Network Node Positions
    grid on;
    axis([-500 500 -500 500])
    title('Sensor Network Geometry'); xlabel('Cross-Range (km)');ylabel('Down-Range (km)')
    [xs,ys] = ginput();
    nSensors = length(xs);
    close;
    pos_sensor = [xs.' ; ys.']*1e3;
    
    
    %% Initialize Emitter Location
    % Define the emitter positions interactively
    grid on
    axis([-500 500 -500 500])
    title('Emitter Position'); xlabel('Cross-Range (km)');ylabel('Down-Range (km)')
    [xe,ye] = ginput();
    pos_emitter = [xe.'; ye.']*1e3;
    close;
    nEmitter = length(xe);
    
    %% Initial estimate of the emitters position
    
%     grid on;
%     axis([-500 500 -500 500])
%     [x_init,y_init] = ginput(nEmitter); % km
%     emitter_init = [x_init.'; y_init.']*1e3;
%     close;
%     for ii = 1: nEmitter
        init = zeros(1,2);
        emitter_init(1,1:2) = init;
%     end
else
    pos_sensor = sensor_pos;
    pos_emitter = emitter_pos; 
    nSensors = size(pos_sensor,2);
    nEmitter = size(pos_emitter,2);
    emitter_init = sensor_pos(:,nSensors);
end

%% Calculate range of each sensor to source and create TOA array for given time
% Initialize array for nEmitters and nSensors
range = zeros(nSensors,nEmitter);
reported_toa = zeros(nSensors,nEmitter);
% reported_toa = cell(1,nEmitter);
for n = 1:nEmitter
    range(:,n) = vecnorm(pos_sensor - pos_emitter(:,n));
    reported_toa(1:nSensors,n) = (range2time(range(:,n),3e8)*0.5).';
end
end

% LOG:
% 1. Date created (updated): 220922 (220928)