function [reported_tdoa, range] = fn_tdoaGen(sensor_pos,emitter_pos,timingErr,tdoaPairs)
%% Generates TOA at sensors of a network for multilateration testing.

% Author: Abeer Nasir Chaudhry®

%% Calculate range of each sensor to source and create TOA array for given time
% Initialize array for nEmitters and nSensors
    c = physconst('Lightspeed');
    nSensor = size(sensor_pos,2);
    nEmitter = size(emitter_pos,2);
    range = zeros(nSensor,nEmitter);
    reported_toa = zeros(nSensor,nEmitter);
    
    for n = 1:nEmitter
        range(:,n) = vecnorm(sensor_pos - emitter_pos(:,n));
        reported_toa(1:nSensor,n) = (range2time(range(:,n),c)*0.5).';
    end
    % Add the measurement noise to the ToAs. The RMS error is the same as
    % std dev when mean of the distribution is 0.
    reported_toa = reported_toa + timingErr*randn(size(reported_toa));
        
    for i = 1:size(tdoaPairs,1)
        reported_tdoa(i) = reported_toa(1) - reported_toa(i+1);
    end
    
end
% LOG:
% 1. Date created (updated): 250110 ()