function [GDOP,axis1,axis2] = gen_GDOPfn(sensor_pos,plane,timingerr,grd_res,axis1_min,axis1_max,axis2_min,axis2_max)
% gdop_fn 

% And computes the Geometric dilution of precision (GDOP) by initializing
% possible source positions over a volume and computing the ratio of error
% in position to the error in range for ach pair of sensors.

% DOP represents the positional precision due to the angular separation
% between the sensors that are used to estimate a source's position.

% Log at end
% Author: Abeer Nasir Chaudhry®

c = 3e8;

tdoaAccuracy = timingerr;

nSensors = size(sensor_pos,2);
sensor_pos(3,1:nSensors) = zeros(1,nSensors);

% Create span
axis1 = linspace(axis1_min,axis1_max,grd_res);
axis2 = linspace(axis2_min,axis2_max,grd_res);

[X,Y] = meshgrid(axis1,axis2);
Z = ones(size(X))*2000;

% All possible source positions
targetPos = [X(:) Y(:) Z(:)]';

%% Decentralised sensor pairing
rxPairs = nchoosek(1:nSensors,2);
n = size(rxPairs,1);

%% Initialize matrices
if strcmp(plane,'xy') == 1 || strcmp(plane,'yx') == 1
    invcova11 = zeros(1,numel(X));
    invcova22 = zeros(1,numel(X));
    invcova12 = zeros(1,numel(X));
    GDOP = zeros(size(X));
else
    invcova11 = zeros(1,numel(Z));
    invcova22 = zeros(1,numel(Z));
    invcova12 = zeros(1,numel(Z));
    GDOP = zeros(size(Z));
end


%% Compute GDOP
for i = 1:n
    % Access a sensor pair from the network
    s1 = sensor_pos(:,rxPairs(i,1));
    s2 = sensor_pos(:,rxPairs(i,2));
    
    % Calculate directional derivative 
    diff_x1 = targetPos(1,:) - s1(1);
    diff_y1 = targetPos(2,:) - s1(2);
    diff_z1 = targetPos(3,:) - s1(3);
    diff_x2 = targetPos(1,:) - s2(1);
    diff_y2 = targetPos(2,:) - s2(2);
    diff_z2 = targetPos(3,:) - s2(3);

    % Range, sensor to target
    r1 = sqrt(diff_x1.^2 + diff_y1.^2 + diff_z1.^2);
    r2 = sqrt(diff_x2.^2 + diff_y2.^2 + diff_z2.^2);
    
    % Directional derivative
    tx = (diff_x1./r1 - diff_x2./r2)./c;
    ty = (diff_y1./r1 - diff_y2./r2)./c;
    % Compute the inverse of Cov mtx for the 
    invcova11 = invcova11 + tx.^2./tdoaAccuracy;
    invcova22 = invcova22 + ty.^2./tdoaAccuracy;
    invcova12 = invcova12 + tx.*ty./tdoaAccuracy;
end

% Get the trace of the covariance mtx
cova11 = invcova22./(invcova11.*invcova22 - invcova12.^2);
cova22 = invcova11./(invcova11.*invcova22 - invcova12.^2);


GDOP(:) = (sqrt(cova11 + cova22));

end


% LOG:
% 1. Date created (updated): 221110 
