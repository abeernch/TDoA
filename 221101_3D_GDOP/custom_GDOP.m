clear;clc

% Sensor Position
% sensor_pos = [00e3 -20e3 15e3 00e3;
%               00e3 15e3 15e3 -15e3;
%               00e3 0.01e3 0.001e3 0.05e3];

sensor_pos = [0	3116.05483775860	-11496.8490106353	7453.21614636302
              0	-7246.75408749039	-21.4179927785758	7273.60074591122
              0  0	                 0	                0];

% No. of sensors
nSensors = size(sensor_pos,2);

% Directional Cosines Mtx intitialize
[hx, hy, hz] = deal(zeros(nSensors,10000));
H = zeros(nSensors,nSensors,10000); 

% Tx speed
c =3e8;

% Sensor timing error
timingError = 1e-9;                 % Sensor Timing Error
rangeError = timingError*c;         % Corresponding range error
% Croa = rangeError^2*eye(nSensors);  % Range Error Covariance matrix

axis1 = linspace(-100e3,100e3,1e2);
axis2 = linspace(-100e3,100e3,1e2);
% axis3 = linspace(-00e3,00e3,1e1);
[X,Y] = meshgrid(axis1,axis2);
Z = zeros(size(X));
% All possible source positions
targetPos = [X(:) Y(:) Z(:)]';

    for i = 1:nSensors
        %% Calculate psudeo ranges to each possible target position from each sensor
        R = vecnorm(targetPos(:,:) - sensor_pos(:,i)) + rangeError;
        %% Compute the directional derivatives
        d_x = sensor_pos(1,i) - targetPos(1,:);
        d_y = sensor_pos(2,i) - targetPos(2,:);
        d_z = sensor_pos(3,i) - targetPos(3,:);

        hx(i,:) = d_x./R;
        hy(i,:) = d_y./R;
        hz(i,:) = d_z./R;

    end
warning('off')
    %% Assign the directional cosines to a 3D matrix
for i = 1:nSensors
    for j = 1:nSensors
        if i == 1
            H(j,i,:) = hx(j,:);
        elseif i == 2
            H(j,i,:) = hy(j,:);            
        elseif i == 3
            H(j,i,:) = hz(j,:);            
        elseif i == 4
            H(j,i,:) = 1; 
        end
    end
end

%% Compute the DOP matrix and various objective DOPs
DOP = pagemtimes(pagetranspose(H),H);
for i = 1:length(DOP)
    DOP(:,:,i) = inv(DOP(:,:,i));
    GDOP(i) = sqrt(sum(diag(DOP(:,:,i))));
    VDOP(i) = DOP(3,3,i);
    HDOP(i) = sqrt((DOP(1,1,i) + DOP(2,2,i) + DOP(3,3,i)));
end

%% Convert the precisions into 2D matrices for vis as contour maps
gdop = reshape(real(GDOP),sqrt(length(GDOP)),sqrt(length(GDOP)));
vdop = reshape(real(VDOP),sqrt(length(GDOP)),sqrt(length(GDOP)));
hdop = reshape(real(HDOP),sqrt(length(GDOP)),sqrt(length(GDOP)));

%% Plotting
figure;[cc,h] = contour(axis1/1e3,axis2/1e3,(gdop)./1e3,[0.1 0.2 0.5 1 2 3 4 5 10 15 20],'LineColor','k');
clabel(cc,h)
hold on
scatter(sensor_pos(1,:)/1e3,sensor_pos(2,:)/1e3,"filled",'MarkerFaceColor','red')
