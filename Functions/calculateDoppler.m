function fd = calculateDoppler(sourcePos,sourceVel,sensorPos,sensorVel,f)
% fd = dop(x1,v1,x2,v2,f)
%
% Given source and sensor at position x1 and x2 with velocity v1 and v2,
% compute the Doppler velocity shift
%
% INPUTS
%   sourcePos      Position vector of N sources (nDim x N), in m
%   sourceVel      Velocity vector of N sources (nDim x N), in m/s
%   sensorPos      Position vector of M sensors (nDim x M), in m
%   sensorVel      Velocity vector of M sensors (nDim x M), in m/s
%   f              Carrier frequency, in Hertz
%
% OUTPUTS
%   fd      Doppler shift for each source, sensor pair (N x M), in Hertz

% Abeer Nasir Chaudhry
% Oct 2024


% % Reshape inputs
[nDim,~] = size(sourcePos);
[nDim2,~] = size(sensorPos);


if nDim~=nDim2
    fprintf('Error: input dimensions do not match.');
    fd = [];
    return
end

% Unit vector from the source to the sensors
u12 = (sensorPos-sourcePos)./sqrt(sum(abs(sensorPos-sourcePos).^2,1));
u21 = -u12; % Unit vector frome sensor to source (will be in the direction 
% opposite to that of the source to sensor)

% Projection of the source velocity on the i-th sensor along the respective
% unit vector and vice versa
vv1 = sum(sourceVel.*u12,1); % source to sensor projection
vv2 = sum(sensorVel.*u21,1); % sensor to source projection. The sensor to 
% source projection is not really needed in our case since we have
% stationary sensors.

% Sum of combined velocity
v = vv1 + vv2;

% Convert to Doppler
c = physconst('LightSpeed');
fd = f .* (1+ v./c);