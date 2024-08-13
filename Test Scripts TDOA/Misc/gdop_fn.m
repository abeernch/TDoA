function [invGDOP,x,y] = gdop_fn(sensor_pos,timingerr)
% gdop_fn 
% The function takes the following inputs:
% 1. Sensor Positions
% 2. Limits for x and y
% 3. TDOA Accuracy (timing error)
% And computes the Geometric dilution of precision (GDOP) by initializing
% possible source positions over a volume and computing the ratio of error
% in position to the error in range for ach pair of sensors.

% DOP represents the positional precision due to the angular separation
% between the sensors that are used to estimate a source's position.

% Log at end
% Author: Abeer Nasir Chaudhry®

c= 3e8;
tdoaAccuracy = timingerr;
nSensors = size(sensor_pos,2);

% Create expanse
x = linspace(-200e3,200e3,1e3);
y = linspace(-200e3,200e3,1e3);
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

% All possible source positions
targetPos = [X(:) Y(:) Z(:)]';

%% Decentralised sensor pairing
rxPairs = nchoosek(1:nSensors,2);
n = size(rxPairs,1);

%% Initialize matrices
invcovxx = zeros(1,numel(X));
invcovyy = zeros(1,numel(X));
invcovxy = zeros(1,numel(X));
invGDOP = zeros(size(X));

%% Compute GDOP
for i = 1:n
    dr1 = sensor_pos(:,rxPairs(i,1));
    dr2 = sensor_pos(:,rxPairs(i,2));
    dx1 = targetPos(1,:) - dr1(1);
    dy1 = targetPos(2,:) - dr1(2);
    dz1 = targetPos(3,:) - dr1(3);
    dx2 = targetPos(1,:) - dr2(1);
    dy2 = targetPos(2,:) - dr2(2);
    dz2 = targetPos(3,:) - dr2(3);

    r1 = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
    r2 = sqrt(dx2.^2 + dy2.^2 + dz2.^2);
    
    tx = (dx1./r1 - dx2./r2)./c;
    ty = (dy1./r1 - dy2./r2)./c;
    
    invcovxx = invcovxx + tx.^2./tdoaAccuracy;
    invcovyy = invcovyy + ty.^2./tdoaAccuracy;
    invcovxy = invcovxy + tx.*ty./tdoaAccuracy;
end

covxx = invcovyy./(invcovxx.*invcovyy - invcovxy.^2);
covyy = invcovxx./(invcovxx.*invcovyy - invcovxy.^2);
invGDOP(:) = (sqrt(covxx + covyy));
end

% LOG:
% 1. Date created (updated): 221004
% 2. Final offset updated to cater for 3D sensor positions